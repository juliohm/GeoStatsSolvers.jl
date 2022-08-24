@testset "IDW" begin
  # basic problem
  geodata = georef((variable=[1.,0.,1.],), [25. 50. 75.;  25. 75. 50.])
  domain  = CartesianGrid(100,100)
  problem = EstimationProblem(geodata, domain, :variable)

  solver = IDW(:variable => (neighbors=3,))

  solution = solve(problem, solver)

  if visualtests
    @test_reference "data/idw.png" contourf(solution,size=(800,400))
  end

  # haversine distance
  geodata = georef((variable=[4.0,-1.0,3.0],), [50.0 100.0 200.0; -30.0 30.0 10.0])
  domain  = CartesianGrid((1.0, -89.0), (359.0, 89.0), dims=(200, 100))
  problem = EstimationProblem(geodata, domain, :variable)

  solver = IDW(:variable => (neighbors=3, distance=Haversine(1.0)))

  solution = solve(problem, solver)

  if visualtests
    @test_reference "data/idw-haversine.png" contourf(solution,size=(800,200))
  end

  # -------------------
  # COMPOSITIONAL DATA
  # -------------------

  table = (z=[Composition(0.1,0.2),Composition(0.3,0.4),Composition(0.5,0.6)],)
  coord = [(25.,25.), (50.,75.), (75.,50.)]
  data = georef(table, coord)
  grid = CartesianGrid(100, 100)

  problem = EstimationProblem(data, grid, :z)

  solver = IDW()

  Random.seed!(2021)
  solution = solve(problem, solver)

  inds = LinearIndices(size(grid))
  S = solution[:z]

  # basic checks
  @test aitchison(S[inds[25,25]], Composition(0.1,0.2)) < 1e-2
  @test aitchison(S[inds[50,75]], Composition(0.3,0.4)) < 1e-2
  @test aitchison(S[inds[75,50]], Composition(0.5,0.6)) < 1e-2
end
