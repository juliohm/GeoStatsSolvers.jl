using GeoStatsSolvers
using Meshes
using GeoTables
using Variography
using GeoStatsBase
using GeoStatsModels
using StatsLearnModels
using Distances
using Distributions
using CategoricalArrays
using Unitful
using CoDa
using LinearAlgebra
using DelimitedFiles
using Test, Random

# dummy definitions
include("dummy.jl")

datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  "ui.jl",
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
