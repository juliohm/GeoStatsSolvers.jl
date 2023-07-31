@testset "LWR" begin
  # 1D regression
  rng = MersenneTwister(2017)
  N = 100
  x = range(0, stop=1, length=N)
  y = x .^ 2 .+ [i / 1000 * randn(rng) for i in 1:N]

  sdata = georef((y=y,), reshape(x, 1, length(x)))
  sdomain = CartesianGrid((0.0,), (1.0,), dims=(N,))
  problem = EstimationProblem(sdata, sdomain, :y)

  solver = LWR(:y => (neighbors=10,))

  sol = solve(problem, solver)

  yhat = sol.y
  yvar = sol.y_variance

  # 2D regression
  sdata = georef((z=[1.0, 0.0, 1.0, 0.0],), [25.0 50.0 75.0 75.0; 25.0 75.0 50.0 25.0])
  sdomain = CartesianGrid(100, 100)
  problem = EstimationProblem(sdata, sdomain, :z)

  solver₃ = LWR(:z => (neighbors=3,))
  solver₄ = LWR(:z => (neighbors=4,))

  sol₃ = solve(problem, solver₃)
  sol₄ = solve(problem, solver₄)

  # Haversine distance
  A = readdlm(joinpath(datadir, "coords.txt"))
  x, y, z = A[:, 1], A[:, 2], A[:, 3]
  sdata = georef((z=z,), [x y]')
  sdomain = begin
    dims = (180, 91)
    start = (1.0, -89.01098901098901)
    finish = (359.0, 89.01098901098901)
    CartesianGrid(start, finish, dims=dims)
  end
  problem = EstimationProblem(sdata, sdomain, :z)

  solver = LWR(:z => (distance=Haversine(6371.0), neighbors=49))

  sol = solve(problem, solver)

  # units
  geodata = georef((; T=[1.0, 0.0, 1.0]u"K"), [25.0 50.0 75.0; 25.0 75.0 50.0])
  domain = CartesianGrid(5, 5)
  problem = EstimationProblem(geodata, domain, :T)
  sol = solve(problem, LWR())
  @test GeoStatsSolvers.elunit(sol.T) == u"K"
  @test GeoStatsSolvers.elunit(sol.T_variance) == u"K^2"

  # affine units
  geodata = georef((; T=[-272.15, -273.15, -272.15]u"°C"), [25.0 50.0 75.0; 25.0 75.0 50.0])
  domain = CartesianGrid(5, 5)
  problem = EstimationProblem(geodata, domain, :T)
  sol = solve(problem, LWR())
  @test GeoStatsSolvers.elunit(sol.T) == u"K"
  @test GeoStatsSolvers.elunit(sol.T_variance) == u"K^2"
end
