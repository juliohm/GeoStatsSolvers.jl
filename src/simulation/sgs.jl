# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SGS(var₁=>param₁, var₂=>param₂, ...)

Sequential Gaussian simulation.

## Parameters

* `variogram` - Variogram model (default to `GaussianVariogram()`)
* `mean`      - Simple Kriging mean
* `degree`    - Universal Kriging degree
* `drifts`    - External Drift Kriging drift functions

Latter options override former options. For example, by specifying
`drifts`, the user is telling the algorithm to ignore `degree` and
`mean`. If no option is specified, Ordinary Kriging is used by
default with the `variogram` only.

* `neighborhood` - Neighborhood on which to search neighbors
* `maxneighbors` - Maximum number of neighbors (default to 10)
* `path`         - Simulation path (default to `LinearPath()`)
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
  @param mean = nothing
  @param degree = nothing
  @param drifts = nothing
  @param neighborhood
  @param minneighbors = 1
  @param maxneighbors = 10
  @param path = nothing
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

      # determine which Kriging variant to use
      if varparams.drifts ≠ nothing
        estimator = ExternalDriftKriging(varparams.variogram, varparams.drifts)
      elseif varparams.degree ≠ nothing
        estimator = UniversalKriging(varparams.variogram, varparams.degree, embeddim(pdomain))
      elseif varparams.mean ≠ nothing
        estimator = SimpleKriging(varparams.variogram, varparams.mean)
      else
        estimator = OrdinaryKriging(varparams.variogram)
      end

      # determine marginal distribution
      marginal = Normal()

      # determine simulation path
      path = varparams.path ≠ nothing ? varparams.path : RandomPath()

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
