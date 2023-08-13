# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Kriging(var₁=>param₁, var₂=>param₂, ...)

The Kriging estimation solver introduced by Matheron 1971.

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
* `path`         - The path algorithm used to iterate over the domain (default to `LinearPath()`)

The `maxneighbors` option can be used to perform approximate Kriging
with a subset of measurements per estimation location. If `maxneighbors`
is not provided, then all measurements are used.

Two `neighborhood` search methods are available:

* If a `neighborhood` is provided, local estimation is performed 
  by sliding the `neighborhood` in the domain.

* If a `neighborhood` is not provided, the estimation is performed 
  using `maxneighbors` nearest neighbors according to `distance`.

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

### References

* Matheron 1971. *The Theory of Regionalized Variables and Its Applications.*
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
  @param path = LinearPath()
end

function preprocess(problem::EstimationProblem, solver::Kriging)
  # retrieve problem info
  pdata = data(problem)
  dtable = values(pdata)
  ddomain = domain(pdata)
  pdomain = domain(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # adjust unit
      cols = Tables.columns(dtable)
      vals = Tables.getcolumn(cols, var)
      z = uadjust(vals)

      # find non-missing samples for variable
      inds = findall(!ismissing, z)

      # assert at least one sample is non-missing
      if isempty(inds)
        throw(AssertionError("all samples of $var are missing, aborting..."))
      end

      # subset of non-missing samples
      vtable = (; var => collect(skipmissing(z)))
      vdomain = view(ddomain, inds)
      samples = georef(vtable, vdomain)

      # determine which Kriging variant to use
      estimator = kriging_ui(pdomain, varparams.variogram, varparams.mean, varparams.degree, varparams.drifts)

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      # determine bounded search method
      searcher = searcher_ui(vdomain, varparams.maxneighbors, varparams.distance, varparams.neighborhood)

      # determine the path algorithm
      path = varparams.path

      # save preprocessed input
      preproc[var] = (; samples, estimator, minneighbors, maxneighbors, searcher, path)
    end
  end

  preproc
end

function solve(problem::EstimationProblem, solver::Kriging)
  # retrieve problem info
  pdomain = domain(problem)

  # preprocess user input
  preproc = preprocess(problem, solver)

  # results for each variable
  μs = []
  σs = []
  for var in name.(variables(problem))
    # maximum number of neighbors
    maxneighbors = preproc[var].maxneighbors

    # non-missing samples
    samples = preproc[var].samples

    # problem for variable
    prob = EstimationProblem(samples, pdomain, var)

    # exact vs. approximate Kriging
    if isnothing(maxneighbors)
      # perform Kriging with all samples as neighbors
      varμ, varσ = exactsolve(prob, var, preproc)
    else
      # perform Kriging with reduced number of neighbors
      varμ, varσ = approxsolve(prob, var, preproc)
    end

    push!(μs, var => varμ)
    push!(σs, Symbol(var, "_variance") => varσ * elunit(varμ)^2)
  end

  georef((; μs..., σs...), pdomain)
end

function exactsolve(problem::EstimationProblem, var::Symbol, preproc)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # unpack preprocessed parameters
  estimator = preproc[var].estimator
  path = preproc[var].path

  # fit estimator once
  krig = fit(estimator, pdata)

  # predict everywhere
  inds = traverse(pdomain, path)
  pred = [predict(krig, var, pdomain[ind]) for ind in inds]

  varμ = first.(pred)
  varσ = last.(pred)

  varμ, varσ
end

function approxsolve(problem::EstimationProblem, var::Symbol, preproc)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # unpack preprocessed parameters
  estimator = preproc[var].estimator
  minneighbors = preproc[var].minneighbors
  maxneighbors = preproc[var].maxneighbors
  searcher = preproc[var].searcher
  path = preproc[var].path

  # pre-allocate memory for neighbors
  neighbors = Vector{Int}(undef, maxneighbors)

  # predict location by location
  inds = traverse(pdomain, path)
  pred = map(inds) do ind
    # centroid of estimation
    center = centroid(pdomain, ind)

    # find neighbors with data
    nneigh = search!(neighbors, center, searcher)

    # skip if there are too few neighbors
    if nneigh < minneighbors
      missing, missing
    else
      # final set of neighbors
      nview = view(neighbors, 1:nneigh)

      # view neighborhood with data
      samples = view(pdata, nview)

      # fit estimator to data
      krig = fit(estimator, samples)

      # save mean and variance
      predict(krig, var, pdomain[ind])
    end
  end

  varμ = first.(pred)
  varσ = last.(pred)

  varμ, varσ
end
