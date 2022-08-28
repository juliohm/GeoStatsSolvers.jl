using GeoStatsSolvers
using Meshes
using GeoStatsBase
using Variography
using Distances
using Distributions
using CategoricalArrays
using CoDa
using LinearAlgebra
using Plots; gr(size=(600,400))
using GeoStatsPlots # TODO: replace by GeoStatsViz
using ReferenceTests, ImageIO
using Test, Random

# learning models from MLJ
using MLJ: @load
dtree = @load DecisionTreeClassifier pkg=DecisionTree verbosity=0

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# dummy definitions
include("dummy.jl")

# list of tests
testfiles = [
  "estimation/idw.jl",
  "estimation/lwr.jl",
  "estimation/kriging.jl",

  "simulation/lu.jl",
  "simulation/fft.jl",
  "simulation/seq.jl",
  "simulation/sgs.jl",
  "simulation/spde.jl",
  "simulation/cookie.jl",

  "learning/pointwise.jl"
]

@testset "GeoStatsSolvers.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
