@testset "CookieCutter" begin
  problem = SimulationProblem(CartesianGrid(100,100), (:facies => Int, :property => Float64), 3)

  solver = CookieCutter(DummySimSolver(:facies=>NamedTuple()), Dict(0=>DummySimSolver(), 1=>DummySimSolver()))

  @test sprint(show, solver) == "CookieCutter"
  @test sprint(show, MIME"text/plain"(), solver) == "CookieCutter\n  └─facies ⇨ DummySimSolver\n    └─0 ⇨ DummySimSolver\n    └─1 ⇨ DummySimSolver"

  Random.seed!(1234)
  solution = solve(problem, solver)

  if visualtests
    @test_reference "data/COOKIE.png" plot(solution)
  end
end
