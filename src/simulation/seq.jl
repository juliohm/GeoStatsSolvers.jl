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
* `mapping`      - Data mapping method (default to `NearestMapping()`)

## Global parameters

* `rng` - random number generator
"""
@simsolver SeqSim begin
  @param estimator
  @param marginal
  @param path = LinearPath()
  @param minneighbors = 1
  @param maxneighbors = 10
  @param neighborhood = nothing
  @param distance = Euclidean()
  @param mapping = NearestMapping()
  @global rng = Random.GLOBAL_RNG
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

      # determine bounded search method
      bsearcher = searcher_ui(pdomain, varparams.maxneighbors, varparams.distance, varparams.neighborhood)

      # determine data mappings
      vmappings = if hasdata(problem)
        map(pdata, pdomain, (var,), varparams.mapping)[var]
      else
        Dict{Int,Int}()
      end

      # save preprocessed input
      preproc[var] = (
        estimator=varparams.estimator,
        minneighbors=minneighbors,
        maxneighbors=maxneighbors,
        marginal=varparams.marginal,
        path=varparams.path,
        bsearcher=bsearcher,
        mappings=vmappings
      )
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
    estimator, minneighbors, maxneighbors, marginal, path, bsearcher, mappings = preproc[var]

    # determine value type
    V = mactypeof[var]

    # pre-allocate memory for result
    realization = Vector{V}(undef, nelements(pdomain))

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # keep track of simulated locations
    simulated = falses(nelements(pdomain))
    if hasdata(problem)
      table = values(pdata)
      cols = Tables.columns(table)
      vals = Tables.getcolumn(cols, var)
      for (ind, datind) in mappings
        realization[ind] = vals[datind]
        simulated[ind] = true
      end
    end

    # simulation loop
    for ind in traverse(pdomain, path)
      if !simulated[ind]
        # coordinates of neighborhood center
        pₒ = centroid(pdomain, ind)

        # find neighbors with previously simulated values
        nneigh = search!(neighbors, pₒ, bsearcher, mask=simulated)

        # choose between marginal and conditional distribution
        if nneigh < minneighbors
          # draw from marginal
          realization[ind] = rand(rng, marginal)
        else
          # final set of neighbors
          ninds = view(neighbors, 1:nneigh)

          # neighborhood with data
          𝒩 = let
            dom = view(pset, ninds)
            val = view(realization, ninds)
            tab = (; var => val)
            georef(tab, dom)
          end

          # fit estimator to data
          fitted = fit(estimator, 𝒩)

          if status(fitted)
            # retrieve element
            uₒ = pdomain[ind]

            # local conditional distribution
            conditional = predictprob(fitted, var, uₒ)

            # draw from conditional
            realization[ind] = rand(rng, conditional)
          else
            # draw from marginal
            realization[ind] = rand(rng, marginal)
          end
        end

        # mark location as simulated and continue
        simulated[ind] = true
      end
    end

    var => realization
  end

  Dict(varreals)
end
