# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    LUGS(var₁=>param₁, var₂=>param₂, ...)

The LU Gaussian simulation solver introduced by Alabert 1987.
The full covariance matrix is built to include all locations
of the simulation domain, and samples from the multivariate
Gaussian are drawn via LU factorization.

## Parameters

* `variogram`     - Theoretical variogram (default to `GaussianVariogram()`)
* `mean`          - Mean of unconditional simulation (default to `0`)
* `factorization` - Factorization method (default to `cholesky`)

## Joint parameters

* `correlation` - correlation coefficient between two covariates (default to `0`).

## Global parameters

* `init` - Data initialization method (default to `NearestInit()`)
* `rng`  - Random number generator (default to `Random.GLOBAL_RNG`)

## Examples

Simulate two variables `var₁` and `var₂` independently:

```julia
julia> LUGS(
         :var₁ => (; variogram=SphericalVariogram(), mean=10.0),
         :var₂ => (; variogram=GaussianVariogram())
       )
```

Simulate two correlated variables `var₁` and `var₂` with correlation `0.7`:

```julia
julia> LUGS(
         :var₁ => (; variogram=SphericalVariogram(), mean=10.0),
         :var₂ => (; variogram=GaussianVariogram()),
         (:var₁, :var₂) => (; correlation=0.7)
       )
```

### References

* Alabert 1987. [The practice of fast conditional simulations
  through the LU decomposition of the covariance matrix]
  (https://link.springer.com/article/10.1007/BF00897191)

* Oliver 2003. [Gaussian cosimulation: modeling of the cross-covariance]
  (https://link.springer.com/article/10.1023%2FB%3AMATG.0000002984.56637.ef)

### Notes

* The solver is only adequate for domains with relatively small
  number of elements (e.g. 100x100 grids) where it is feasible to
  factorize the full covariance.

* For larger domains (e.g. 3D grids), other solvers are preferred
  such as [`SGS`](@ref) and [`FFTGS`](@ref).
"""
@simsolver LUGS begin
  @param variogram = GaussianVariogram()
  @param mean = nothing
  @param factorization = cholesky
  @jparam correlation = 0.0
  @global init = NearestInit()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::LUGS)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)
  pvars = variables(problem)

  # retrieve global parameters
  init = solver.init

  # initialize buffers for realizations and simulation mask
  buff, mask = initbuff(pdomain, pvars, init, data=pdata)

  # result of preprocessing
  preproc = Dict()

  for covars in covariables(problem, solver)
    conames = covars.names
    coparams = Dict()

    # 1 or 2 variables can be simulated simultaneously
    @assert length(conames) ∈ (1, 2) "invalid number of covariables"

    # preprocess parameters for individual variables
    for var in conames
      # get parameters for variable
      varparams = covars.params[Set([var])]

      # determine variogram model
      γ = varparams.variogram

      # determine factorization method
      fact = varparams.factorization

      # check stationarity
      @assert isstationary(γ) "variogram model must be stationary"

      # retrieve data locations and data values in domain
      dlocs = findall(mask[var])
      z₁ = view(buff[var], dlocs)

      # retrieve simulation locations
      slocs = [l for l in 1:nelements(pdomain) if l ∉ dlocs]

      # create views of the domain
      𝒟d = [centroid(pdomain, i) for i in dlocs]
      𝒟s = [centroid(pdomain, i) for i in slocs]

      # covariance between simulation locations
      C₂₂ = sill(γ) .- Variography.pairwise(γ, 𝒟s)

      if isempty(dlocs)
        d₂ = zero(eltype(z₁))
        L₂₂ = fact(Symmetric(C₂₂)).L
      else
        # covariance beween data locations
        C₁₁ = sill(γ) .- Variography.pairwise(γ, 𝒟d)
        C₁₂ = sill(γ) .- Variography.pairwise(γ, 𝒟d, 𝒟s)

        L₁₁ = fact(Symmetric(C₁₁)).L
        B₁₂ = L₁₁ \ C₁₂
        A₂₁ = B₁₂'

        d₂ = A₂₁ * (L₁₁ \ z₁)
        L₂₂ = fact(Symmetric(C₂₂ - A₂₁ * B₁₂)).L
      end

      if !isnothing(varparams.mean) && !isempty(dlocs)
        @warn "mean can only be specified in unconditional simulation"
      end

      # mean for unconditional simulation
      μ = isnothing(varparams.mean) ? zero(eltype(z₁)) : varparams.mean

      # save preprocessed parameters for variable
      coparams[Set([var])] = (z₁, d₂, L₂₂, μ, dlocs, slocs)
    end

    # preprocess joint parameters
    if length(conames) == 2
      # get parameters for pair of variables
      jparams = covars.params[conames]

      # 0-lag correlation between variables
      ρ = jparams.correlation

      # save preprocessed parameters for pair of variables
      coparams[conames] = ρ
    end

    push!(preproc, conames => coparams)
  end

  preproc
end

function solvesingle(::SimulationProblem, covars::NamedTuple, solver::LUGS, preproc)
  # random number generator
  rng = solver.rng

  # preprocessed parameters
  conames = covars.names
  coparams = preproc[conames]

  # plain list of variable names
  vars = collect(conames)

  # simulate first variable
  v₁ = first(vars)
  Y₁, w₁ = lusim(rng, coparams[Set([v₁])])
  result = Dict(v₁ => Y₁)

  # simulate second variable
  if length(conames) == 2
    ρ = coparams[conames]
    v₂ = last(vars)
    Y₂, w₂ = lusim(rng, coparams[Set([v₂])], ρ, w₁)
    push!(result, v₂ => Y₂)
  end

  result
end

function lusim(rng, params, ρ=nothing, w₁=nothing)
  # unpack parameters
  z₁, d₂, L₂₂, μ, dlocs, slocs = params

  # number of points in domain
  npts = length(dlocs) + length(slocs)

  # allocate memory for result
  y = Vector{eltype(z₁)}(undef, npts)

  # conditional simulation
  w₂ = randn(rng, size(L₂₂, 2))
  if isnothing(ρ)
    y₂ = d₂ .+ L₂₂ * w₂
  else
    y₂ = d₂ .+ L₂₂ * (ρ * w₁ + √(1 - ρ^2) * w₂)
  end

  # hard data and simulated values
  y[dlocs] = z₁
  y[slocs] = y₂

  # adjust mean in case of unconditional simulation
  isempty(dlocs) && (y .+= μ)

  y, w₂
end
