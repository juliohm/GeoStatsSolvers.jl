# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SGS(var₁=>param₁, var₂=>param₂, ...)

Sequential Gaussian simulation.

## Parameters

* `variogram`    - Variogram model (default to `GaussianVariogram()`)
* `mean`         - mean for simple Kriging (default to `0.0`)
* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `mapping`      - Data mapping method (default to `NearestMapping()`)

For each location in the simulation `path`, a maximum number of
neighbors `maxneighbors` is used to fit a Gaussian distribution.
The neighbors are searched according to a `neighborhood`.

## Global parameters

* `rng` - random number generator (default to `Random.GLOBAL_RNG`)

### References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)
"""
@simsolver SGS begin
  @param variogram = GaussianVariogram()
  @param mean = 0.0
  @param path = LinearPath()
  @param minneighbors = 1
  @param maxneighbors = 10
  @param neighborhood = nothing
  @param distance = Euclidean()
  @param mapping = NearestMapping()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::SGS)
  # retrieve problem info
  pdomain = domain(problem)

  params = []
  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine simple Kriging estimator
      estimator = SimpleKriging(varparams.variogram, varparams.mean)

      # determine marginal distribution
      μ = varparams.mean
      σ = √sill(varparams.variogram)
      marginal = Normal(μ, σ)

      # determine simulation path
      path = isnothing(varparams.path) ? RandomPath() : varparams.path

      # determine data mapping
      mapping = varparams.mapping

      # equivalent parameters for SeqSim solver
      param = (estimator=estimator,
               neighborhood=varparams.neighborhood,
               minneighbors=varparams.minneighbors,
               maxneighbors=varparams.maxneighbors,
               marginal=marginal, path=path,
               mapping=mapping)

      push!(params, var => param)
    end
  end

  preprocess(problem, SeqSim(params...))
end

solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::SGS, preproc) =
  solvesingle(problem, covars, SeqSim(rng=solver.rng), preproc)
