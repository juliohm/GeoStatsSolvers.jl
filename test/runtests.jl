using GeoStatsSolvers
using Meshes
using GeoTables
using GeoStatsBase
using Variography
using GeoStatsModels
using StatsLearnModels
using Distances
using Distributions
using CategoricalArrays
using LinearAlgebra
using DelimitedFiles
using Unitful
using Random
using CoDa
using Test

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
