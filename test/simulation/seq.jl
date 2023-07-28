@testset "SeqSim" begin
  sdata = georef((z=[1.,0.,1.],), [(25.,25.), (50.,75.), (75.,50.)])
  sgrid = CartesianGrid(100,100)

  prob1 = SimulationProblem(sgrid, :z => Float64, 3)
  prob2 = SimulationProblem(sdata, sgrid, :z, 3)

  rng = MersenneTwister(123)
  solver = SeqSim(:z => (estimator=DummyEstimator(), marginal=Normal()),
                  rng=rng)

  usol = solve(prob1, solver)
  csol = solve(prob2, solver)
end
