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
  dtable = values(pdata)
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
      dcols = Tables.columns(dtable)
      dvals = Tables.getcolumn(dcols, var)
      dinds = findall(!ismissing, dvals)
      𝒮 = view(pdata, dinds)
      𝒟 = domain(𝒮)
      𝒯 = values(𝒮)
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
      X = [coordinates(centroid(𝒟, i)) for i in 1:n]
      if D isa NearestNeighbors.MinkowskiMetric
        tree = KDTree(X, D)
      else
        tree = BallTree(X, D)
      end

      # adjust unit
      cols = Tables.columns(𝒯)
      vals = Tables.getcolumn(cols, var)
      unit = elunit(vals)
      z = uadjust(unit, vals)

      # estimation loop
      inds = traverse(pdomain, LinearPath())
      pred = map(inds) do ind
        x = coordinates(centroid(pdomain, ind))

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

        ẑₒ, r̂ₒ
      end

      varμ = first.(pred)
      varσ = last.(pred)

      push!(μs, var => urevert(unit, varμ))
      push!(σs, Symbol(var, "_variance") => varσ * absoluteunit(unit)^2)
    end
  end

  georef((; μs..., σs...), pdomain)
end
