import GeoStatsBase: fit, predict, predictprob, status

# ----------------
# DUMMY ESTIMATOR
# ----------------

struct DummyEstimator end
struct FittedDummyEstimator end

fit(::DummyEstimator, data) = FittedDummyEstimator()
predict(::FittedDummyEstimator, var, pₒ) = (0, 1)
predictprob(::FittedDummyEstimator, var, pₒ) = Normal(0, 1)
status(::FittedDummyEstimator) = true

# ------------------------
# DUMMY SIMULATION SOLVER
# ------------------------

import GeoStatsBase: solvesingle

@simsolver DummySimSolver begin end
function solvesingle(problem::SimulationProblem, covars::NamedTuple, ::DummySimSolver, preproc)
  npts = nelements(domain(problem))
  reals = map(collect(covars.names)) do var
    V = variables(problem)[var]
    real = vcat(fill(zero(V), npts ÷ 2), fill(one(V), npts ÷ 2))
    var => real
  end
  Dict(reals)
end
