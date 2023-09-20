@testset "Kriging" begin
  # -----------------
  # 1D PROBLEM (ALL)
  # -----------------

  data1D = georef((; z=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]), collect(0.0:10.0:100.0)')
  grid1D = CartesianGrid(100)
  problem = EstimationProblem(data1D, grid1D, :z)

  global_kriging = KrigingSolver(:z => (; variogram=GaussianVariogram(range=35.0, nugget=0.0)))
  nearest_kriging = KrigingSolver(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3))
  local_kriging = KrigingSolver(
    :z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3, neighborhood=MetricBall(100.0))
  )

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

  solver = KrigingSolver(:z => (; variogram=GaussianVariogram(range=35.0, nugget=0.0)))

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

  solver = KrigingSolver(:z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3))

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

  solver = KrigingSolver(
    :z => (variogram=GaussianVariogram(range=35.0, nugget=0.0), maxneighbors=3, neighborhood=MetricBall(100.0))
  )

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

  solver = KrigingSolver(
    :z => (
      variogram=GaussianVariogram(range=35.0, nugget=0.0),
      maxneighbors=3,
      neighborhood=MetricBall(100.0),
      path=MultiGridPath()
    )
  )

  Random.seed!(2021)
  sol = solve(problem, solver)
end
