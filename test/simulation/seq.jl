@testset "SeqSim" begin
  sdata = georef((z=[1.0, 0.0, 1.0],), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  sgrid = CartesianGrid(100, 100)

  prob1 = SimulationProblem(sgrid, :z => Float64, 3)
  prob2 = SimulationProblem(sdata, sgrid, :z, 3)

  rng = MersenneTwister(123)
  solver = SeqSim(:z => (estimator=DummyEstimator(), marginal=Normal()), rng=rng)

  usol = solve(prob1, solver)
  csol = solve(prob2, solver)
end
