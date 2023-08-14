# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SeqSim(var₁=>param₁, var₂=>param₂, ...)

A sequential simulation solver.

For each location in the simulation `path`, a maximum number
of neighbors `maxneighbors` is used to fit a distribution.
The neighbors are searched according to a `neighborhood`,
and in case there are none, use a `marginal` distribution.

## Parameters

* `estimator`    - CDF estimator
* `marginal`     - Marginal distribution
* `path`         - Simulation path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

## Global parameters

* `init` - Data initialization method (default to `NearestInit()`)
* `rng`  - Random number generator (default to `Random.GLOBAL_RNG`)
"""
@simsolver SeqSim begin
  @param estimator
  @param marginal
  @param path = LinearPath()
  @param minneighbors = 1
  @param maxneighbors = 10
  @param neighborhood = nothing
  @param distance = Euclidean()
  @global init = NearestInit()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::SeqSim)
  # retrieve problem info
  pdomain = domain(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[Set([var])]

      # extract paramaters
      estimator = varparams.estimator
      marginal = varparams.marginal
      path = varparams.path
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors
      neighborhood = varparams.neighborhood
      distance = varparams.distance

      # determine search method
      searcher = searcher_ui(pdomain, maxneighbors, distance, neighborhood)

      # save preprocessed input
      preproc[var] = (
        estimator,
        marginal,
        path,
        minneighbors,
        maxneighbors,
        searcher=searcher
      )
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::SeqSim, preproc)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)
  pvars = variables(problem)

  # retrieve global parameters
  init = solver.init
  rng = solver.rng

  # initialize buffers for realization and simulation mask
  buff, mask = initbuff(pdomain, pvars, init, data=pdata)

  # consider point set with centroids for now
  pset = PointSet([centroid(pdomain, ind) for ind in 1:nelements(pdomain)])

  varreals = map(collect(covars.names)) do var
    # unpack preprocessed parameters
    estimator, marginal, path, minneighbors, maxneighbors, searcher = preproc[var]

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # retrieve realization and mask for variable
    realization = buff[var]
    simulated = mask[var]

    # simulation loop
    for ind in traverse(pdomain, path)
      if !simulated[ind]
        # search neighbors with simulated data
        nneigh = search!(neighbors, pset[ind], searcher, mask=simulated)

        if nneigh < minneighbors
          # draw from marginal
          realization[ind] = rand(rng, marginal)
        else
          # neighborhood with data
          neigh = let
            ijk = view(neighbors, 1:nneigh)
            dom = view(pset, ijk)
            val = view(realization, ijk)
            tab = (; var => val)
            georef(tab, dom)
          end

          # fit distribution estimator
          fitted = fit(estimator, neigh)

          # draw from conditional or marginal
          distribution = if status(fitted)
            predictprob(fitted, var, pset[ind])
          else
            marginal
          end
          realization[ind] = rand(rng, distribution)
        end

        # mark location as simulated and continue
        simulated[ind] = true
      end
    end

    var => realization
  end

  Dict(varreals)
end
