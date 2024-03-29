# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CookieCutter(master, others)

A cookie-cutter simulation solver.

## Parameters

* `master` - Master simulation solver (a.k.a. facies solver)
* `others` - Solvers for each value of the master realization

## Examples

Simulate lithofacies with image quilting and fill porosity
values with LU Gaussian simulation:

```julia
julia> f  = IQ(:facies => (TI=IMG, template=(30,30,1)))
julia> p₀ = LUGS(:poro => (variogram=SphericalVariogram(range=10.),))
julia> p₁ = LUGS(:poro => (variogram=SphericalVariogram(range=20.),))
julia> CookieCutter(f, Dict(0=>p₀, 1=>p₁))
```
"""
struct CookieCutter <: SimulationSolver
  master::SimulationSolver
  others::Dict{Number,SimulationSolver}
end

function solve(problem::SimulationProblem, solver::CookieCutter)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)
  preals = nreals(problem)

  # master variable
  pvars = variables(problem)
  mvars = targets(solver.master)
  @assert length(mvars) == 1 "one single variable must be specified in master solver"
  mvar = mvars[1]
  @assert mvar ∈ keys(pvars) "invalid variable in master solver"
  mtype = pvars[mvar]

  # other variables
  ovars = Tuple(var => V for (var, V) in pairs(pvars) if var ≠ mvar)
  @assert length(ovars) > 0 "cookie-cutter requires problem with more than one target variable"

  # define master and others problem
  if !isnothing(pdata)
    mproblem = SimulationProblem(pdata, pdomain, mvar, preals)
    oproblem = SimulationProblem(pdata, pdomain, first.(ovars), preals)
  else
    mproblem = SimulationProblem(pdomain, mvar => mtype, preals)
    oproblem = SimulationProblem(pdomain, ovars, preals)
  end

  # realizations of master variable
  msol = solve(mproblem, solver.master)
  mreals = msol[mvar]

  # pre-allocate memory for realizations
  reals = Dict(var => [Vector{V}(undef, nelements(pdomain)) for i in 1:preals] for (var, V) in ovars)

  # solve other problems
  for (mval, osolver) in solver.others
    osol = solve(oproblem, osolver)

    # cookie-cutter step
    for (var, V) in ovars
      vreals = osol[var]
      for i in 1:preals
        mask = mreals[i] .== mval
        reals[var][i][mask] .= vreals[i][mask]
      end
    end
  end

  realizations = merge(reals, Dict(mvar => mreals))

  Ensemble(pdomain, realizations)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, solver::CookieCutter)
  print(io, "CookieCutter")
end

function Base.show(io::IO, ::MIME"text/plain", solver::CookieCutter)
  mvar = targets(solver.master)[1]
  println(io, solver)
  println(io, "  └─", mvar, " ⇨ ", solver.master)
  lines = ["    └─$val ⇨ $osolver" for (val, osolver) in solver.others]
  print(io, join(lines, "\n"))
end
