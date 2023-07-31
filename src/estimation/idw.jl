# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    IDW(var₁=>param₁, var₂=>param₂, ...)

The inverse distance weighting estimation solver introduced
in the very early days of geostatistics by Shepard 1968. It
is perhaps the simplest first attempt in the literature to
perform estimation based on the notion of proximity to data
locations.

## Parameters

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `nothing`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance defined in Distances.jl (default to `Euclidean()`)
* `power`        - Power of the distances (default to `1`)

The `maxneighbors` option can be used to perform IDW estimation
with a subset of data points per estimation location. If `maxneighbors`
is set to `nothing` then all data points will be used. Two neighborhood
search methods are available depending on the value of `neighborhood`:

* If a `neighborhood` is provided, local estimation is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is provided, the estimation is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

### References

Shepard 1968. *A two-dimensional interpolation function for irregularly-spaced data.*

### Notes

* This implementation makes use of k-d trees from the NearestNeighbors.jl package,
  which leads to a fast estimation method for spatial domains with large number of
  elements.

* Although this method is recommended for fast assessment of a new field, it has
  poor statistical properties (lacks covariance model) and should mainly be used
  for qualitative purposes.
"""
@estimsolver IDW begin
  @param minneighbors = 1
  @param maxneighbors = nothing
  @param neighborhood = nothing
  @param distance = Euclidean()
  @param power = 1
end

function solve(problem::EstimationProblem, solver::IDW)
  # retrieve problem info
  pdata = data(problem)
  dtable = values(pdata)
  pdomain = domain(problem)

  # result for each variable
  μs = []
  σs = []

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # retrieve non-missing data
      dcols = Tables.columns(dtable)
      dvals = Tables.getcolumn(dcols, var)
      dinds = findall(!ismissing, dvals)
      𝒮 = view(pdata, dinds)
      𝒟 = domain(𝒮)
      𝒯 = values(𝒮)
      n = nelements(𝒟)

      # determine power of distances
      p = varparams.power

      # determine distance
      distance = varparams.distance

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      @assert n > 0 "estimation requires data"

      @assert p > 0 "power must be positive"

      @assert maxneighbors ≤ n "invalid number of maxneighbors"

      @assert minneighbors < maxneighbors "invalid number of minneighbors"

      # determine bounded search method
      bsearcher = searcher_ui(𝒟, maxneighbors, distance, varparams.neighborhood)

      # pre-allocate memory for neighbors
      neighbors = Vector{Int}(undef, maxneighbors)

      # pre-compute the centroid coordinates
      X = [coordinates(centroid(𝒟, i)) for i in 1:n]

      # adjust unit
      cols = Tables.columns(𝒯)
      vals = Tables.getcolumn(cols, var)
      z = uadjust(vals)

      # estimation loop
      inds = traverse(pdomain, LinearPath())
      pred = map(inds) do ind
        # centroid of estimation
        center = centroid(pdomain, ind)

        # find neighbors with data
        nneigh = search!(neighbors, center, bsearcher)

        # skip if there are too few neighbors
        if nneigh < minneighbors
          missing, missing
        else
          x = coordinates(center)
          is = view(neighbors, 1:nneigh)
          ds = [distance(x, X[i]) for i in is]
          ws = 1 ./ ds .^ p
          Σw = sum(ws)

          if isinf(Σw) # some distance is zero?
            j = findfirst(iszero, ds)
            μ = z[is[j]]
            σ = zero(eltype(ds))
          else
            ws /= Σw
            vs = view(z, is)
            μ = sum(ws[i] * vs[i] for i in eachindex(vs))
            σ = minimum(ds)
          end

          μ, σ
        end
      end

      varμ = first.(pred)
      varσ = last.(pred)

      push!(μs, var => varμ)
      push!(σs, Symbol(var, "_distance") => varσ)
    end
  end

  georef((; μs..., σs...), pdomain)
end
