@testset "FFTGS" begin
  # isotropic simulation
  Random.seed!(2019)
  problem = SimulationProblem(CartesianGrid(100,100), :z=>Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  sol     = solve(problem, solver)

  if visualtests
    @test_reference "data/FFT-iso.png" plot(sol,size=(900,300))
  end

  # anisotropic simulation
  Random.seed!(2019)
  problem = SimulationProblem(CartesianGrid(100,100), :z=>Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(MetricBall((20.,5.))),))
  sol = solve(problem, solver)

  if visualtests
    @test_reference "data/FFT-aniso.png" plot(sol,size=(900,300))
  end

  # simulation on view of grid
  Random.seed!(2022)
  grid    = CartesianGrid(100,100)
  vgrid   = view(grid, 1:5000)
  problem = SimulationProblem(vgrid, :z=>Float64, 3)
  solver  = FFTGS(:z => (variogram=GaussianVariogram(range=10.),))
  sol = solve(problem, solver)
  @test domain(sol) == vgrid
  @test length(sol[1].z) == 5000
end
