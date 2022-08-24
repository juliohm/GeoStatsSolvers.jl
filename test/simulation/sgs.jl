@testset "SGS" begin
  𝒮 = georef((z=[1.,0.,1.],), [25. 50. 75.; 25. 75. 50.])
  𝒟 = CartesianGrid((100,100), (0.5,0.5), (1.0,1.0))
  N = 3

  problem₁ = SimulationProblem(𝒮, 𝒟, :z, N)
  problem₂ = SimulationProblem(𝒟, :z=>Float64, N)

  solver = SGS(:z => (variogram=SphericalVariogram(range=35.),
                      neighborhood=MetricBall(30.)))

  Random.seed!(2017)
  sol₁ = solve(problem₁, solver)
  sol₂ = solve(problem₂, solver)

  # basic checks
  reals = sol₁[:z]
  inds = LinearIndices(size(𝒟))
  @test all(reals[i][inds[25,25]] == 1. for i in 1:N)
  @test all(reals[i][inds[50,75]] == 0. for i in 1:N)
  @test all(reals[i][inds[75,50]] == 1. for i in 1:N)
end
