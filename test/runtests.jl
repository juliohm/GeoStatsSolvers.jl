using GeoStatsSolvers
using Meshes
using GeoStatsBase
using Variography
using Distances
using Distributions
using CategoricalArrays
using Unitful
using CoDa
using MLJ: @load
using LinearAlgebra
using DelimitedFiles
using Test, Random

# dummy definitions
include("dummy.jl")

datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  "estimation/idw.jl",
  "estimation/lwr.jl",
  "estimation/krig.jl",

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
