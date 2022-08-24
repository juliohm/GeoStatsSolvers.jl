# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Kriging(var₁=>param₁, var₂=>param₂, ...)

A polyalgorithm Kriging estimation solver.

Each pair `var=>param` specifies the `KrigingParam` `param`
for the Kriging variable `var`. In order to avoid boilerplate
code, the constructor expects pairs of `Symbol` and `NamedTuple`
instead.

## Parameters

* `variogram` - Variogram model (default to `GaussianVariogram()`)
* `mean`      - Simple Kriging mean
* `degree`    - Universal Kriging degree
* `drifts`    - External Drift Kriging drift functions

Latter options override former options. For example, by specifying
`drifts`, the user is telling the algorithm to ignore `degree` and
`mean`. If no option is specified, Ordinary Kriging is used by
default with the `variogram` only.

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `nothing`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

The `maxneighbors` option can be used to perform approximate Kriging
with a subset of data points per estimation location. If `maxneighbors`
is set to `nothing`, global Kriging is performed. Two neighborhood
search methods are available depending on the value of `neighborhood`:

* If a `neighborhood` is provided, local Kriging is performed by sliding
  the `neighborhood` in the domain.

* If `neighborhood` is not provided, the Kriging system is built using
  `maxneighbors` nearest neighbors according to a `distance`.

## Examples

Solve the variable `:var₁` with Simple Kriging by specifying
the `mean`, and the variable `:var₂` with Universal Kriging
by specifying the `degree` and the `variogram` model.

```julia
julia> Kriging(
  :var₁ => (mean=1.,),
  :var₂ => (degree=1, variogram=SphericalVariogram(range=20.))
)
```

Solve all variables of the problem with the default parameters
(i.e. Ordinary Kriging with unit Gaussian variogram):

```julia
julia> Kriging()
```
"""
@estimsolver Kriging begin
  @param variogram = GaussianVariogram()
  @param mean = nothing
  @param degree = nothing
  @param drifts = nothing
  @param minneighbors = 1
  @param maxneighbors = nothing
  @param neighborhood = nothing
  @param distance = Euclidean()
end

function preprocess(problem::EstimationProblem, solver::Kriging)
  # retrieve problem info
  pdata = data(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine which Kriging variant to use
      estimator = kriging_ui(domain(pdata),
                             varparams.variogram,
                             varparams.mean,
                             varparams.degree,
                             varparams.drifts)

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      # determine bounded search method
      bsearcher = searcher_ui(domain(pdata),
                              varparams.maxneighbors,
                              varparams.distance,
                              varparams.neighborhood)

      # save preprocessed input
      preproc[var] = (estimator=estimator,
                      minneighbors=minneighbors,
                      maxneighbors=maxneighbors,
                      bsearcher=bsearcher)
    end
  end

  preproc
end

function solve(problem::EstimationProblem, solver::Kriging)
  # preprocess user input
  preproc = preprocess(problem, solver)

  # results for each variable
  μs = []; σs = []
  for var in name.(variables(problem))
    if preproc[var].maxneighbors ≠ nothing
      # perform Kriging with reduced number of neighbors
      varμ, varσ = solve_approx(problem, var, preproc)
    else
      # perform Kriging with all data points as neighbors
      varμ, varσ = solve_exact(problem, var, preproc)
    end

    push!(μs, var => varμ)
    push!(σs, Symbol(var,"_variance") => varσ)
  end

  georef((; μs..., σs...), domain(problem))
end

function solve_approx(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher = preproc[var]

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # predict at all locations
    locations = traverse(pdomain, LinearPath())
    predictions = map(locations) do loc

      pₒ = centroid(pdomain, loc)

      # find neighbors with previously estimated values
      nneigh = search!(neighbors, pₒ, bsearcher)

      # skip location in there are too few neighbors
      if nneigh < minneighbors
        missing, missing
      else
        # final set of neighbors
        nview = view(neighbors, 1:nneigh)

        # view neighborhood with data
        𝒟 = view(pdata, nview)

        # fit estimator to data
        krig = fit(estimator, 𝒟)

        # retrieve element at location
        uₒ = pdomain[loc]

        # save mean and variance
        predict(krig, var, uₒ)
      end
    end

    varμ = first.(predictions)
    varσ = last.(predictions)

    varμ, varσ
end

function solve_exact(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher = preproc[var]

    # retrieve non-missing data
    locs = findall(!ismissing, pdata[var])
    𝒟 = view(pdata, locs)

    # fit estimator once
    krig = fit(estimator, 𝒟)

    # predict at all locations
    locations = traverse(pdomain, LinearPath())
    predictions = map(locations) do loc
      uₒ = pdomain[loc]
      predict(krig, var, uₒ)
    end

    varμ = first.(predictions)
    varσ = last.(predictions)

    varμ, varσ
end
