@testset "Kriging" begin
  # -----------------
  # 1D PROBLEM (ALL)
  # -----------------

  data1D = georef((z=[0.0,0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1,0.0],), collect(0.0:10.0:100.0)')
  grid1D = CartesianGrid(100)
  problem = EstimationProblem(data1D, grid1D, :z)

  global_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),))
  nearest_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),
                                   maxneighbors=3))
  local_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),
                                 maxneighbors=3, neighborhood=MetricBall(100.)))

  Random.seed!(2021)
  solvers   = [global_kriging, nearest_kriging, local_kriging]
  solnames  = ["global", "nearest", "local"]
  solutions = [solve(problem, solver) for solver in solvers]

  if visualtests
    for i in 1:3
      solution, sname = solutions[i], solnames[i]
      @test_reference "data/krig-1D-$(sname).png" plot(solution, size=(600,200))
    end
  end

  # --------------------
  # 2D PROBLEM (GLOBAL)
  # --------------------

  data2D = georef((z=[1.,0.,1.],), [(25.,25.), (50.,75.), (75.,50.)])
  grid2D = CartesianGrid((100,100), (0.5,0.5), (1.0,1.0))

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),))

  Random.seed!(2021)
  solution = solve(problem, solver)
  Z = asarray(solution, :z)
  @test isapprox(Z[25,25], 1., atol=1e-3)
  @test isapprox(Z[50,75], 0., atol=1e-3)
  @test isapprox(Z[75,50], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-global.png" contourf(solution, size=(800,400))
  end

  # ---------------------
  # 2D PROBLEM (NEAREST)
  # ---------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),))

  Random.seed!(2021)
  solution = solve(problem, solver)
  Z = asarray(solution, :z)
  @test isapprox(Z[25,25], 1., atol=1e-3)
  @test isapprox(Z[50,75], 0., atol=1e-3)
  @test isapprox(Z[75,50], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-nearest.png" contourf(solution, size=(800,400))
  end

  # -------------------
  # 2D PROBLEM (LOCAL)
  # -------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),
                          maxneighbors=3, neighborhood=MetricBall(100.)))

  Random.seed!(2021)
  solution = solve(problem, solver)

  # basic checks
  inds = LinearIndices(size(grid2D))
  S = solution[:z]
  @test isapprox(S[inds[25,25]], 1., atol=1e-3)
  @test isapprox(S[inds[50,75]], 0., atol=1e-3)
  @test isapprox(S[inds[75,50]], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-local.png" contourf(solution, size=(800,400))
  end

  # -------------------
  # COMPOSITIONAL DATA
  # -------------------

  table = (z=[Composition(0.1,0.2),Composition(0.3,0.4),Composition(0.5,0.6)],)
  coord = [(25.,25.), (50.,75.), (75.,50.)]
  data = georef(table, coord)
  grid = CartesianGrid(100, 100)

  problem = EstimationProblem(data, grid, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.),))

  Random.seed!(2021)
  solution = solve(problem, solver)

  inds = LinearIndices(size(grid))
  S = solution[:z]

  # basic checks
  @test aitchison(S[inds[25,25]], Composition(0.1,0.2)) < 1e-3
  @test aitchison(S[inds[50,75]], Composition(0.3,0.4)) < 1e-3
  @test aitchison(S[inds[75,50]], Composition(0.5,0.6)) < 1e-3
end
