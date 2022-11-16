# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    FFTGS(var₁=>param₁, var₂=>param₂, ...)

FFT Gaussian simulation.

## Parameters

* `variogram` - theoretical variogram (default to `GaussianVariogram()`)
* `mean`      - mean of Gaussian field (default to `0`)

In the case of conditional simulation, the following parameters
can be passed to the underlying [`Kriging`](@ref) solver:

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

## Global parameters

* `threads` - number of threads in FFT (default to all physical cores)
* `rng` - random number generator (default to `Random.GLOBAL_RNG`)

### References

* Gutjahr 1997. [General joint conditional simulations using a fast
  Fourier transform method](https://link.springer.com/article/10.1007/BF02769641)

* Gómez-Hernández, J. & Srivastava, R. 2021. [One Step at a Time: The Origins
  of Sequential Simulation and Beyond](https://link.springer.com/article/10.1007/s11004-021-09926-0)
"""
@simsolver FFTGS begin
  @param  variogram = GaussianVariogram()
  @param  mean = 0
  @param  minneighbors = 1
  @param  maxneighbors = 10
  @param  neighborhood = nothing
  @param  distance = Euclidean()
  @global threads = cpucores()
  @global rng = Random.GLOBAL_RNG
end

function preprocess(problem::SimulationProblem, solver::FFTGS)
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
      if !isstationary(γ)
        throw(ArgumentError("variogram model must be stationary"))
      end

      # compute covariances between centroid and all points
      𝒟c = [centroid(pgrid, cindex)]
      𝒟p = [centroid(pgrid, eindex) for eindex in 1:nelms]
      cs = sill(γ) .- pairwise(γ, 𝒟c, 𝒟p)
      C  = reshape(cs, dims)

      # move to frequency domain
      F = sqrt.(abs.(fft(fftshift(C))))
      F[1] = zero(V) # set reference level

      # perform Kriging in case of conditional simulation
      z̄, krig, dinds = nothing, nothing, nothing
      if hasdata(problem)
        pdata   = data(problem)
        dtable  = values(pdata)
        ddomain = domain(pdata)
        if var ∈ Tables.schema(dtable).names
          # estimate conditional mean
          kdat = georef(dtable, centroid.(ddomain))
          kdom = PointSet(centroid.(pdomain))
          prob = EstimationProblem(kdat, kdom, var)
          krig = Kriging(var => (
              variogram = γ, mean = μ,
              minneighbors = varparams.minneighbors,
              maxneighbors = varparams.maxneighbors,
              neighborhood = varparams.neighborhood,
              distance     = varparams.distance
          ))
          ksol = solve(prob, krig)
          z̄ = getproperty(ksol, var)

          # find data locations in problem domain
          ndata    = nelements(ddomain)
          point(i) = centroid(ddomain, i)
          searcher = KNearestSearch(pdomain, 1)
          found    = [search(point(i), searcher) for i in 1:ndata]
          dinds    = unique(first.(found))
        end
      end

      # save preprocessed inputs for variable
      preproc[var] = (γ=γ, μ=μ, F=F, z̄=z̄, krig=krig, dinds=dinds)
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
    γ, μ, F, z̄, krig, dinds = preproc[var]

    # determine value type
    V = mactypeof[var]

    # perturbation in frequency domain
    P = F .* exp.(im .* angle.(fft(rand(rng, V, dims))))

    # move back to time domain
    Z = real(ifft(P))

    # adjust mean and variance
    σ² = Statistics.var(Z, mean=zero(V))
    Z .= √(sill(γ) / σ²) .* Z .+ μ

    # unconditional realization
    zᵤ = Z[inds]

    # perform conditioning if necessary
    z = if isnothing(krig)
      zᵤ # we are all set
    else
      # view realization at data locations
      dtable  = (;var => view(zᵤ, dinds))
      ddomain = view(pdomain, dinds)

      # solve estimation problem
      kdat = georef(dtable, centroid.(ddomain))
      kdom = PointSet(centroid.(pdomain))
      prob = EstimationProblem(kdat, kdom, var)
      ksol = solve(prob, krig)
      z̄ᵤ   = getproperty(ksol, var)

      # add residual field
      z̄ .+ (zᵤ .- z̄ᵤ)
    end

    var => z
  end

  Dict(varreal)
end
