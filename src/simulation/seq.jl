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
* `neighborhood` - Geospatial neighborhood
* `maxneighbors` - Maximum number of neighbors
* `marginal`     - Marginal distribution
* `path`         - Simulation path
* `mapping`      - Data mapping method

## Global parameters

* `rng` - random number generator
"""
@simsolver SeqSim begin
  @param estimator
  @param neighborhood
  @param minneighbors
  @param maxneighbors
  @param marginal
  @param path
  @param mapping
  @global rng
end

function preprocess(problem::SimulationProblem, solver::SeqSim)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      # determine neighbor search method
      neigh     = varparams.neighborhood
      searcher  = BallSearch(pdomain, neigh)
      bsearcher = BoundedSearch(searcher, maxneighbors)

      # determine data mappings
      vmappings = if hasdata(problem)
        map(pdata, pdomain, (var,), varparams.mapping)[var]
      else
        Dict()
      end

      # save preprocessed input
      preproc[var] = (estimator=varparams.estimator,
                      minneighbors=minneighbors,
                      maxneighbors=maxneighbors,
                      marginal=varparams.marginal,
                      path=varparams.path,
                      bsearcher=bsearcher,
                      mappings=vmappings)
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::SeqSim, preproc)
  # random number generator
  rng = solver.rng

  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # compute variogram between centroids
  pset = PointSet(centroid.(pdomain))

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  varreals = map(covars.names) do var
    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors,
    marginal, path, bsearcher, mappings = preproc[var]

    # determine value type
    V = mactypeof[var]

    # pre-allocate memory for result
    realization = Vector{V}(undef, nelements(pdomain))

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # keep track of simulated locations
    simulated = falses(nelements(pdomain))
    for (loc, datloc) in mappings
      realization[loc] = pdata[var][datloc]
      simulated[loc] = true
    end

    # simulation loop
    for location in traverse(pdomain, path)
      if !simulated[location]
        # coordinates of neighborhood center
        pₒ = centroid(pdomain, location)

        # find neighbors with previously simulated values
        nneigh = search!(neighbors, pₒ, bsearcher, mask=simulated)

        # choose between marginal and conditional distribution
        if nneigh < minneighbors
          # draw from marginal
          realization[location] = rand(rng, marginal)
        else
          # final set of neighbors
          nview = view(neighbors, 1:nneigh)

          # view neighborhood with data
          tab = (; var => view(realization, nview))
          dom = view(pset, nview)
          𝒟   = georef(tab, dom)

          # fit estimator to data
          fitted = fit(estimator, 𝒟)

          if status(fitted)
            # retrieve element
            uₒ = pdomain[location]

            # local conditional distribution
            conditional = predictprob(fitted, var, uₒ)

            # draw from conditional
            realization[location] = rand(rng, conditional)
          else
            # draw from marginal
            realization[location] = rand(rng, marginal)
          end
        end

        # mark location as simulated and continue
        simulated[location] = true
      end
    end

    var => realization
  end

  Dict(varreals)
end
