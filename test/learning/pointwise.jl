@testset "PointwiseLearn" begin
  rng = MersenneTwister(1234)
  tree = @load DecisionTreeClassifier pkg=DecisionTree verbosity=0

  # synthetic data
  f(x,y) = sin(4*(abs(x)+abs(y))) < 0 ? 1 : 0 
  X = [sin(i/10) for i in 1:100, j in 1:100]
  Y = [sin(j/10) for i in 1:100, j in 1:100]
  Z = categorical(f.(X,Y))
  ϵ₁ = 0.1randn(rng, Float64, size(X))
  ϵ₂ = 0.1randn(rng, Float64, size(Y))

  # source and target data
  S = georef((X=X,Y=Y,Z=Z))
  T = georef((X=X+ϵ₁,Y=Y+ϵ₂))

  # view versions
  inds = shuffle(rng, 1:nitems(S))
  Sv = view(S, inds)
  Tv = view(T, inds)

  # classification task
  𝓉 = ClassificationTask((:X,:Y), :Z)

  # learning problems
  𝒫₁ = LearningProblem(S, T, 𝓉)
  𝒫₂ = LearningProblem(Sv, Tv, 𝓉)

  # pointwise solver
  ℒ = PointwiseLearn(tree())

  R₁ = solve(𝒫₁, ℒ)
  R₂ = solve(𝒫₂, ℒ)

  # error is small
  @test mean(S.Z .!= R₁.Z) < 0.15
  @test mean(Sv.Z .!= R₂.Z) < 0.15
end
