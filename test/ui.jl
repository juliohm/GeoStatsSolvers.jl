@testset "UI" begin
  # ------------
  # searcher UI
  # ------------

  domain = PointSet(rand(2, 3))
  method = GeoStatsSolvers.searcher_ui(domain, 2, Euclidean(), nothing)
  @test method isa KNearestSearch
  @test maxneighbors(method) == 2
  method = GeoStatsSolvers.searcher_ui(domain, 2, nothing, MetricBall(1.0))
  @test method isa KBallSearch
  @test maxneighbors(method) == 2
  method = GeoStatsSolvers.searcher_ui(domain, nothing, Euclidean(), nothing)
  @test method isa KNearestSearch
  @test maxneighbors(method) == 3
  method = @test_logs (:warn, "Invalid maximum number of neighbors. Adjusting to 3...") GeoStatsSolvers.searcher_ui(
    domain,
    4,
    Euclidean(),
    nothing
  )
  @test method isa KNearestSearch
  @test maxneighbors(method) == 3

  # -----------
  # Kriging UI
  # -----------

  grid = CartesianGrid(10, 10)
  krig = GeoStatsSolvers.kriging_ui(grid, GaussianVariogram(), nothing, nothing, nothing)
  @test krig isa OrdinaryKriging
  krig = GeoStatsSolvers.kriging_ui(grid, GaussianVariogram(), 0.0, nothing, nothing)
  @test krig isa SimpleKriging
  krig = GeoStatsSolvers.kriging_ui(grid, GaussianVariogram(), nothing, 2, nothing)
  @test krig isa UniversalKriging
  krig = GeoStatsSolvers.kriging_ui(grid, GaussianVariogram(), nothing, nothing, [x -> 1])
  @test krig isa ExternalDriftKriging
end
