# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SGS(var₁=>param₁, var₂=>param₂, ...)

The sequential Gaussian simulation solver introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using [`Kriging`](@ref), and assigns a value to the center of the
neighborhood by sampling from this distribution.

## Parameters

* `variogram`    - Variogram model (default to `GaussianVariogram()`)
* `mean`         - mean for simple Kriging (default to `0.0`)
* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

For each location in the simulation `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

## Global parameters

* `init` - Data initialization method (default to `NearestInit()`)
* `rng`  - Random number generator (default to `Random.GLOBAL_RNG`)

### References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This solver is very sensitive to the simulation path and number of
  samples. Care must be taken to make sure that neighborhoods have
  enough samples for the Kriging estimator.
"""
@simsolver SGS begin
  @param variogram = GaussianVariogram()
  @param mean = 0.0
  @param path = LinearPath()
  @param minneighbors = 1
  @param maxneighbors = 10
  @param neighborhood = nothing
  @param distance = Euclidean()
  @global init = NearestInit()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::SGS)
  params = []
  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[Set([var])]

      # determine Kriging estimator
      estimator = GeoStatsModels.SimpleKriging(varparams.variogram, varparams.mean)

      # determine marginal distribution
      μ = varparams.mean
      σ = √sill(varparams.variogram)
      marginal = Normal(μ, σ)

      # forward remaining parameters
      path = varparams.path
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors
      neighborhood = varparams.neighborhood
      distance = varparams.distance

      # equivalent parameters for SeqSim solver
      param = (; estimator, marginal, path, minneighbors, maxneighbors, neighborhood, distance)

      push!(params, var => param)
    end
  end

  preprocess(problem, SeqSim(params...))
end

solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::SGS, preproc) =
  solvesingle(problem, covars, SeqSim(init=solver.init, rng=solver.rng), preproc)
