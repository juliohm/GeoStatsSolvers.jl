# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    LWR(var₁=>param₁, var₂=>param₂, ...)

Locally weighted regression estimation solver.

## Parameters

* `neighbors` - Number of neighbors (default to all the data)
* `distance`  - A distance from Distances.jl (default to `Euclidean()`)
* `weightfun` - Weighting function (default to `exp(-3*h^2/2)`)

### References

* Stone 1977. *Consistent non-parametric regression.*
* Cleveland 1979. *Robust locally weighted regression and smoothing scatterplots.*
* Cleveland & Grosse 1991. *Computational methods for local regression.*
"""
@estimsolver LWR begin
  @param neighbors = nothing
  @param distance = Euclidean()
  @param weightfun = h -> exp(-3*h^2)
end

function solve(problem::EstimationProblem, solver::LWR)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  # result for each variable
  μs = []; σs = []

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine value type
      V = mactypeof[var]

      # retrieve non-missing data
      locs = findall(!ismissing, pdata[var])
      𝒟 = view(pdata, locs)
      n = nelements(𝒟)

      # determine number of nearest neighbors to use
      k = isnothing(varparams.neighbors) ? n : varparams.neighbors

      # determine distance type
      D = varparams.distance

      # weight function
      w = varparams.weightfun

      @assert n > 0 "estimation requires data"

      @assert k ≤ n "invalid number of neighbors"

      # fit search tree
      X = [coordinates(centroid(𝒟, i)) for i in 1:nelements(𝒟)]
      if D isa NearestNeighbors.MinkowskiMetric
        tree = KDTree(X, D)
      else
        tree = BallTree(X, D)
      end

      # lookup non-missing values
      z = 𝒟[var]

      # pre-allocate memory for results
      varμ = Vector{V}(undef, nelements(pdomain))
      varσ = Vector{V}(undef, nelements(pdomain))

      # estimation loop
      for loc in traverse(pdomain, LinearPath())
        x = coordinates(centroid(pdomain, loc))

        # find neighbors
        is, ds = knn(tree, x, k)
        δs = ds ./ maximum(ds)

        # weighted least-squares
        Wₗ = Diagonal(w.(δs))
        Xₗ = [ones(eltype(x), k) reduce(hcat, X[is])']
        zₗ = view(z, is)
        θₗ = Xₗ'*Wₗ*Xₗ \ Xₗ'*Wₗ*zₗ

        # linear combination of response values
        xₒ = [one(eltype(x)); x]
        ẑₒ = θₗ ⋅ xₒ
        rₗ = Wₗ*Xₗ*(Xₗ'*Wₗ*Xₗ\xₒ)
        r̂ₒ = norm(rₗ)

        varμ[loc] = ẑₒ
        varσ[loc] = r̂ₒ
      end

      push!(μs, var => varμ)
      push!(σs, Symbol(var,"_variance") => varσ)
    end
  end

  georef((; μs..., σs...), pdomain)
end
