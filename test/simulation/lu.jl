@testset "LUGS" begin
  𝒮 = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
  𝒟 = CartesianGrid(100)

  # ----------------------
  # conditional simulation
  # ----------------------
  problem = SimulationProblem(𝒮, 𝒟, :z, 2)

  rng = MersenneTwister(123)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),), rng=rng)

  solution = solve(problem, solver)

  # ------------------------
  # unconditional simulation
  # ------------------------
  problem = SimulationProblem(𝒟, :z=>Float64, 2)

  rng = MersenneTwister(123)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),), rng=rng)

  solution = solve(problem, solver)

  # -------------
  # co-simulation
  # -------------
  𝒟 = CartesianGrid(500)
  problem = SimulationProblem(𝒟, (:z=>Float64,:y=>Float64), 1)

  rng = MersenneTwister(123)
  solver = LUGS(:z => (variogram=SphericalVariogram(range=10.),),
                :y => (variogram=GaussianVariogram(range=10.),),
                (:z,:y) => (correlation=0.95,),
                rng=rng)

  solution = solve(problem, solver)

  # -----------
  # 2D example
  # -----------
  𝒟 = CartesianGrid(100,100)
  problem = SimulationProblem(𝒟, :z=>Float64, 3)

  rng = MersenneTwister(123)
  solver = LUGS(:z => (variogram=GaussianVariogram(range=10.),), rng=rng)

  solution = solve(problem, solver)

  # -------------------
  # anisotropy example
  # -------------------
  𝒟 = CartesianGrid(100,100)
  problem = SimulationProblem(𝒟, :z=>Float64, 3)

  rng = MersenneTwister(123)
  ball = MetricBall((20.,5.))
  solver = LUGS(:z => (variogram=GaussianVariogram(ball),), rng=rng)

  solution = solve(problem, solver)

  # ---------------------
  # custom factorization
  # ---------------------
  𝒮 = georef((z=[0.,1.,0.,1.,0.],), [0. 25. 50. 75. 100.])
  𝒟 = CartesianGrid(100)
  problem = SimulationProblem(𝒮, 𝒟, :z, 1)

  rng = MersenneTwister(123)
  solver1 = LUGS(:z => (variogram=SphericalVariogram(range=10.),factorization=lu), rng=rng)
  solver2 = LUGS(:z => (variogram=SphericalVariogram(range=10.),factorization=cholesky), rng=rng)

  solution1 = solve(problem, solver1)
  solution2 = solve(problem, solver2)
end
