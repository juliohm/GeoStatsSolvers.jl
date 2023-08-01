@testset "Kriging" begin
  # -----------------
  # 1D PROBLEM (ALL)
  # -----------------

  data1D = georef((; z=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]), collect(0.0:10.0:100.0)')
  grid1D = CartesianGrid(100)
  problem = EstimationProblem(data1D, grid1D, :z)

  global_kriging = Kriging(:z => (; variogram=GaussianVariogram(range=35.0, nugget=0.0)))
  nearest_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3))
  local_kriging =
    Kriging(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3, neighborhood=MetricBall(100.0)))

  Random.seed!(2021)
  solvers = [global_kriging, nearest_kriging, local_kriging]
  solnames = ["global", "nearest", "local"]
  solutions = [solve(problem, solver) for solver in solvers]

  # --------------------
  # 2D PROBLEM (GLOBAL)
  # --------------------

  data2D = georef((; z=[1.0, 0.0, 1.0]), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  grid2D = CartesianGrid((100, 100), (0.5, 0.5), (1.0, 1.0))

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (; variogram=GaussianVariogram(range=35.0, nugget=0.0)))

  Random.seed!(2021)
  sol = solve(problem, solver)
  Z = asarray(sol, :z)
  @test isapprox(Z[25, 25], 1.0, atol=1e-3)
  @test isapprox(Z[50, 75], 0.0, atol=1e-3)
  @test isapprox(Z[75, 50], 1.0, atol=1e-3)

  # ---------------------
  # 2D PROBLEM (NEAREST)
  # ---------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3))

  Random.seed!(2021)
  sol = solve(problem, solver)
  Z = asarray(sol, :z)
  @test isapprox(Z[25, 25], 1.0, atol=1e-3)
  @test isapprox(Z[50, 75], 0.0, atol=1e-3)
  @test isapprox(Z[75, 50], 1.0, atol=1e-3)

  # -------------------
  # 2D PROBLEM (LOCAL)
  # -------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver =
    Kriging(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3, neighborhood=MetricBall(100.0)))

  Random.seed!(2021)
  sol = solve(problem, solver)

  # basic checks
  S = sol.z
  inds = LinearIndices(size(grid2D))
  @test isapprox(S[inds[25, 25]], 1.0, atol=1e-3)
  @test isapprox(S[inds[50, 75]], 0.0, atol=1e-3)
  @test isapprox(S[inds[75, 50]], 1.0, atol=1e-3)

  # -------------------------
  # 2D PROBLEM (CUSTOM PATH)
  # -------------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(
    :z => (
      variogram=GaussianVariogram(range=35.0, nugget=0.0),
      maxneighbors=3,
      neighborhood=MetricBall(100.0),
      path=MultiGridPath()
    )
  )

  Random.seed!(2021)
  sol = solve(problem, solver)

  # units
  geodata = georef((; T=[1.0, 0.0, 1.0]u"K"), [25.0 50.0 75.0; 25.0 75.0 50.0])
  domain = CartesianGrid(5, 5)
  problem = EstimationProblem(geodata, domain, :T)
  sol = solve(problem, Kriging())
  @test GeoStatsSolvers.elunit(sol.T) == u"K"
  @test GeoStatsSolvers.elunit(sol.T_variance) == u"K^2"

  # affine units
  geodata = georef((; T=[-272.15, -273.15, -272.15]u"°C"), [25.0 50.0 75.0; 25.0 75.0 50.0])
  domain = CartesianGrid(5, 5)
  problem = EstimationProblem(geodata, domain, :T)
  sol = solve(problem, Kriging())
  @test GeoStatsSolvers.elunit(sol.T) == u"K"
  @test GeoStatsSolvers.elunit(sol.T_variance) == u"K^2"

  # -------------------
  # COMPOSITIONAL DATA
  # -------------------

  table = (; z=[Composition(0.1, 0.2), Composition(0.3, 0.4), Composition(0.5, 0.6)])
  coord = [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)]
  data = georef(table, coord)
  grid = CartesianGrid(100, 100)

  problem = EstimationProblem(data, grid, :z)

  solver = Kriging(:z => (; variogram=GaussianVariogram(range=35.0)))

  Random.seed!(2021)
  sol = solve(problem, solver)
  inds = LinearIndices(size(grid))
  S = sol.z
  @test aitchison(S[inds[25, 25]], Composition(0.1, 0.2)) < 1e-3
  @test aitchison(S[inds[50, 75]], Composition(0.3, 0.4)) < 1e-3
  @test aitchison(S[inds[75, 50]], Composition(0.5, 0.6)) < 1e-3
end
