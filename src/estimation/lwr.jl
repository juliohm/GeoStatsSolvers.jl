# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    LWRSolver(var₁=>param₁, var₂=>param₂, ...)

The locally weighted regression (a.k.a. LOESS) estimation solver
introduced by Cleveland 1979. It is the most natural generalization
of [`IDWSolver`](@ref) in which one is allowed to use a custom weight
function instead of distance-based weights.

## Parameters

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `nothing`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - A distance from Distances.jl (default to `Euclidean()`)
* `weightfun`    - Weighting function (default to `exp(-3 * h^2)`)
* `path`         - The path algorithm used to iterate over the domain (default to `LinearPath()`)

The `maxneighbors` option can be used to perform locally weighted regression
with a subset of measurements per estimation location. If `maxneighbors`
is not provided, then all measurements are used.

Two `neighborhood` search methods are available:

* If a `neighborhood` is provided, local estimation is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is not provided, the estimation is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

### References

* Stone 1977. *Consistent non-parametric regression.*
* Cleveland 1979. *Robust locally weighted regression and smoothing scatterplots.*
* Cleveland & Grosse 1991. *Computational methods for local regression.*

### Notes

* This implementation makes use of k-d trees from the NearestNeighbors.jl
  package for fast data lookup.

* Locally weighted regression (LWR or LOESS) is a popular non-parametric
  solver, however it still has poor statistical properties compared to
  other estimation solvers such as [`KrigingSolver`](@ref) that explicitly
  model spatial correlation.

* In the current implementation, the estimation variance is computed
  assuming Gaussian residuals. 
"""
@estimsolver LWRSolver begin
  @param minneighbors = 1
  @param maxneighbors = nothing
  @param neighborhood = nothing
  @param distance = Euclidean()
  @param weightfun = h -> exp(-3 * h^2)
  @param path = LinearPath()
end

function solve(problem::EstimationProblem, solver::LWRSolver)
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
      varparams = covars.params[Set([var])]

      # retrieve non-missing data
      dcols = Tables.columns(dtable)
      dvals = Tables.getcolumn(dcols, var)
      dinds = findall(!ismissing, dvals)
      𝒮 = view(pdata, dinds)
      𝒟 = domain(𝒮)
      𝒯 = values(𝒮)
      n = nelements(𝒟)

      # retrieve solver params
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors
      neighborhood = varparams.neighborhood
      distance = varparams.distance
      weightfun = varparams.weightfun
      path = varparams.path

      nmin = minneighbors
      nmax = isnothing(maxneighbors) ? n : min(maxneighbors, n)

      @assert n > 0 "estimation requires data"
      @assert nmin ≤ nmax "invalid min/max number of neighbors"

      # determine bounded search method
      searcher = searcher_ui(𝒟, maxneighbors, distance, neighborhood)

      # pre-allocate memory for neighbors
      neighbors = Vector{Int}(undef, nmax)

      # pre-allocate memory for distances
      distances = Vector{coordtype(𝒟)}(undef, nmax)

      # adjust unit
      cols = Tables.columns(𝒯)
      vals = Tables.getcolumn(cols, var)
      z = uadjust(vals)

      # estimation loop
      inds = traverse(pdomain, path)
      pred = map(inds) do ind
        # centroid of estimation
        center = centroid(pdomain, ind)

        # find neighbors with data
        nneigh = searchdists!(neighbors, distances, center, searcher)

        # skip if there are too few neighbors
        if nneigh < nmin
          missing, missing
        else
          x = coordinates(center)
          is = view(neighbors, 1:nneigh)
          ds = view(distances, 1:nneigh)
          δs = ds ./ maximum(ds)

          # weighted least-squares
          X = mapreduce(i -> coordinates(centroid(𝒟, i)), hcat, is)
          Wₗ = Diagonal(weightfun.(δs))
          Xₗ = [ones(eltype(x), nneigh) X']
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
