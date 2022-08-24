# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SPDEGS(var₁=>param₁, var₂=>param₂, ...)

SPDE Gaussian simulation.

## Parameters

* `sill`  - Sill or total variance (default to `1.0`)
* `range` - Range or correlation length (default to `1.0`)

### References

* Lindgren et al. 2011. [An explicit link between Gaussian fields and
  Gaussian Markov random fields: the stochastic partial differential
  equation approach](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2011.00777.x)
"""
@simsolver SPDEGS begin
  @param sill = 1.0
  @param range = 1.0
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::SPDEGS)
  hasdata(problem) && @error "conditional simulation is not implemented"

  # retrieve problem info
  𝒟 = domain(problem)
  d = paramdim(𝒟)

  # Beltrami-Laplace discretization
  B = laplacematrix(𝒟)
  M = measurematrix(𝒟)
  Δ = inv(M) * B

  # result of preprocessing
  preproc = Dict()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine sill and range
      σ = varparams.sill
      𝓁 = varparams.range

      @assert σ > 0 "sill must be positive"
      @assert 𝓁 > 0 "range must be positive"

      # LHS of SPDE (κ² - Δ)Z = τW
      α = 2one(σ+𝓁)
      ν = α - d / 2
      κ = 1 / 𝓁
      A = κ^2*I - Δ

      # covariance structure
      τ² = σ^2 * κ^(2ν) * (4π)^(d/2) * gamma(α) / gamma(ν)
      Q  = A'A / τ²

      # factorization
      F = cholesky(Array(Q))
      L = inv(Array(F.U))

      # save preprocessed inputs for variable
      preproc[var] = (L=L,)
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::SPDEGS, preproc)
  # random number generator
  rng = solver.rng

  # retrieve problem info
  𝒟 = domain(problem)
  n = nvertices(𝒟)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  # target variables
  vars = covars.names

  # simulation at vertices
  varreal = map(vars) do var
    # unpack preprocessed parameters
    L = preproc[var].L

    # determine value type
    V = mactypeof[var]

    # perform simulation
    w = randn(rng, V, n)
    z = L*w

    var => z
  end

  # vertex table
  vtable = (; varreal...)

  # change of support
  vdata = meshdata(𝒟, vtable=vtable)
  edata = integrate(vdata, vars...)

  # columns of element table
  cols = Tables.columns(values(edata))

  Dict(vars .=> cols)
end