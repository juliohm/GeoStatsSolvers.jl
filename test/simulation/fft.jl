@testset "FFTGS" begin
  # isotropic simulation
  Random.seed!(2019)
  problem = SimulationProblem(CartesianGrid(100,100), :z => Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  sol     = solve(problem, solver)

  # anisotropic simulation
  Random.seed!(2019)
  problem = SimulationProblem(CartesianGrid(100,100), :z => Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(MetricBall((20.,5.))),))
  sol     = solve(problem, solver)

  # simulation on view of grid
  Random.seed!(2022)
  grid    = CartesianGrid(100,100)
  vgrid   = view(grid, 1:5000)
  problem = SimulationProblem(vgrid, :z => Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  sol     = solve(problem, solver)
  @test domain(sol) == vgrid
  @test length(sol[1].z) == 5000

  # conditional simulation
  Random.seed!(2022)
  table   = (z=[1.,-1.,1.],)
  coords  = [(25.,25.), (50.,75.), (75.,50.)]
  samples = georef(table, coords)
  sdomain = CartesianGrid(100,100)
  problem = SimulationProblem(samples, sdomain, :z => Float64, 100)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  sol     = solve(problem, solver)
end
