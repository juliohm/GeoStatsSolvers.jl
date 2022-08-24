@testset "SeqSim" begin
  Random.seed!(1234)
  sdata = georef((z=[1.,0.,1.],), [(25.,25.), (50.,75.), (75.,50.)])
  sgrid = CartesianGrid(100,100)

  prob1 = SimulationProblem(sgrid, :z => Float64, 3)
  prob2 = SimulationProblem(sdata, sgrid, :z, 3)

  rng = MersenneTwister(123)
  solver = SeqSim(:z => (estimator=DummyEstimator(),
                         neighborhood=MetricBall(10.),
                         minneighbors=1, maxneighbors=10,
                         marginal=Normal(), path=RandomPath(),
                         mapping=NearestMapping()), rng=rng)

  usol = solve(prob1, solver)
  csol = solve(prob2, solver)

  # reals = csol[:z]
  # inds = LinearIndices(size(sgrid))
  # @test all(reals[i][inds[25,25]] == 1. for i in 1:3)
  # @test all(reals[i][inds[50,75]] == 0. for i in 1:3)
  # @test all(reals[i][inds[75,50]] == 1. for i in 1:3)
end
