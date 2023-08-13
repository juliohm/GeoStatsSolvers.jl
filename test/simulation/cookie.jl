@testset "CookieCutter" begin
  problem = SimulationProblem(CartesianGrid(100, 100), (:f => Int, :z => Float64), 3)

  solver = CookieCutter(DummySimSolver(:f => NamedTuple()), Dict(0 => DummySimSolver(), 1 => DummySimSolver()))

  @test sprint(show, solver) == "CookieCutter"
  @test sprint(show, MIME"text/plain"(), solver) ==
        "CookieCutter\n  └─f ⇨ DummySimSolver\n    └─0 ⇨ DummySimSolver\n    └─1 ⇨ DummySimSolver"

  Random.seed!(1234)
  solution = solve(problem, solver)
end
