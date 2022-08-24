# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    FFTGS(var₁=>param₁, var₂=>param₂, ...)

FFT Gaussian simulation.

## Parameters

* `variogram` - theoretical variogram (default to `GaussianVariogram()`)
* `mean`      - mean of Gaussian field (default to `0`)

## Global parameters

* `threads` - number of threads in FFT (default to all physical cores)
* `rng` - random number generator (default to `Random.GLOBAL_RNG`)

### References

* Gutjahr 1997. [General joint conditional simulations using a fast
Fourier transform method](https://link.springer.com/article/10.1007/BF02769641)
"""
@simsolver FFTGS begin
  @param  variogram = GaussianVariogram()
  @param  mean = 0.0
  @global threads = cpucores()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::FFTGS)
  hasdata(problem) && @error "conditional simulation is not implemented"
  
  # retrieve problem info
  pdomain  = domain(problem)
  pgrid, _ = unview(pdomain)
  dims     = size(pgrid)
  nelms    = nelements(pgrid)
  center   = CartesianIndex(dims .÷ 2)
  cindex   = LinearIndices(dims)[center]

  # number of threads in FFTW
  FFTW.set_num_threads(solver.threads)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  # result of preprocessing
  preproc = Dict()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine value type
      V = mactypeof[var]

      # determine variogram model and mean
      γ = varparams.variogram
      μ = varparams.mean

      # check stationarity
      @assert isstationary(γ) "variogram model must be stationary"

      # compute covariances between centroid and all points
      𝒟c = [centroid(pgrid, cindex)]
      𝒟p = [centroid(pgrid, eindex) for eindex in 1:nelms]
      covs = sill(γ) .- pairwise(γ, 𝒟c, 𝒟p)
      C = reshape(covs, dims)

      # move to frequency domain
      F = sqrt.(abs.(fft(fftshift(C))))
      F[1] = zero(V) # set reference level

      # save preprocessed inputs for variable
      preproc[var] = (γ=γ, μ=μ, F=F)
    end
  end

  preproc
end

function solvesingle(problem::SimulationProblem, covars::NamedTuple, solver::FFTGS, preproc)
  # random number generator
  rng = solver.rng

  # retrieve problem info
  pdomain     = domain(problem)
  pgrid, inds = unview(pdomain)
  dims        = size(pgrid)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  varreal = map(covars.names) do var
    # unpack preprocessed parameters
    γ, μ, F = preproc[var]

    # determine value type
    V = mactypeof[var]

    # perturbation in frequency domain
    P = F .* exp.(im .* angle.(fft(rand(rng, V, dims))))

    # move back to time domain
    Z = real(ifft(P))

    # adjust mean and variance
    σ² = Statistics.var(Z, mean=zero(V))
    Z .= √(sill(γ) / σ²) .* Z .+ μ

    # flatten result
    var => Z[inds]
  end

  Dict(varreal)
end
