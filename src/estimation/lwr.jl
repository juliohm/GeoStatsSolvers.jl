# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    LWR(var₁=>param₁, var₂=>param₂, ...)

The locally weighted regression (a.k.a. LOESS) estimation solver
introduced by Cleveland 1979. It is the most natural generalization
of [`IDW`](@ref) in which one is allowed to use a custom weight
function instead of distance-based weights.

## Parameters

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `nothing`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance from Distances.jl (default to `Euclidean()`)
* `weightfun`    - Weighting function (default to `exp(-3 * h^2)`)

The `maxneighbors` option can be used to perform IDW estimation
with a subset of data points per estimation location. If `maxneighbors`
is set to `nothing` then all data points will be used. Two neighborhood
search methods are available depending on the value of `neighborhood`:

* If a `neighborhood` is provided, local estimation is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is provided, the estimation is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

### References

* Stone 1977. *Consistent non-parametric regression.*
* Cleveland 1979. *Robust locally weighted regression and smoothing scatterplots.*
* Cleveland & Grosse 1991. *Computational methods for local regression.*

### Notes

* This implementation makes use of k-d trees from the NearestNeighbors.jl
  package for fast data lookup.

* Locally weighted regression (LWR or LOESS) is a popular non-parametric
  method, however it still has poor statistical properties compared to
  other estimation methods such as [`Kriging`](@ref) that explicitly
  model spatial correlation.

* In the current implementation, the estimation variance is computed
  assuming Gaussian residuals. 
"""
@estimsolver LWR begin
  @param minneighbors = 1
  @param maxneighbors = nothing
  @param neighborhood = nothing
  @param distance = Euclidean()
  @param weightfun = h -> exp(-3 * h^2)
end

function solve(problem::EstimationProblem, solver::LWR)
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

      # weight function
      w = varparams.weightfun

      # determine distance
      distance = varparams.distance

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      @assert n > 0 "estimation requires data"
      @assert minneighbors < n "invalid number of minneighbors"
      if !isnothing(maxneighbors)
        @assert maxneighbors ≤ n "invalid number of maxneighbors"
        @assert minneighbors < maxneighbors "invalid number of minneighbors"
      end

      # determine bounded search method
      bsearcher = searcher_ui(𝒟, maxneighbors, distance, varparams.neighborhood)

      # pre-allocate memory for neighbors
      neighbors = Vector{Int}(undef, isnothing(maxneighbors) ? n : maxneighbors)

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
          δs = ds ./ maximum(ds)

          # weighted least-squares
          Wₗ = Diagonal(w.(δs))
          Xₗ = [ones(eltype(x), nneigh) reduce(hcat, X[is])']
          zₗ = view(z, is)
          θₗ = Xₗ' * Wₗ * Xₗ \ Xₗ' * Wₗ * zₗ

          # linear combination of response values
          xₒ = [one(eltype(x)); x]
          ẑₒ = θₗ ⋅ xₒ
          rₗ = Wₗ * Xₗ * (Xₗ' * Wₗ * Xₗ \ xₒ)
          r̂ₒ = norm(rₗ)

          ẑₒ, r̂ₒ
        end
      end

      varμ = first.(pred)
      varσ = last.(pred)

      push!(μs, var => varμ)
      push!(σs, Symbol(var, "_variance") => varσ * elunit(varμ)^2)
    end
  end

  georef((; μs..., σs...), pdomain)
end
