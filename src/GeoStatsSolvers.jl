# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoStatsSolvers

using Meshes
using GeoTables
using Variography
using GeoStatsBase
using GeoStatsModels
using StatsLearnModels

using Tables
using Distributions
using TableTransforms
using Distances: Euclidean
using Bessels: gamma
using Unitful: Units, AffineUnits
using Unitful
using CpuId
using FFTW

using LinearAlgebra
using Statistics
using Random

import GeoStatsBase: preprocess, solve, solvesingle

include("utils.jl")
include("ui.jl")

include("estimation/idw.jl")
include("estimation/lwr.jl")
include("estimation/krig.jl")

include("simulation/lu.jl")
include("simulation/fft.jl")
include("simulation/seq.jl")
include("simulation/sgs.jl")
include("simulation/spde.jl")
include("simulation/cookie.jl")

include("learning/pointwise.jl")

export
  # -----------
  # ESTIMATION
  # -----------
  IDWSolver,
  LWRSolver,
  KrigingSolver,

  # -----------
  # SIMULATION
  # -----------
  # generic solvers
  SeqSim,

  # concrete solvers
  LUGS,
  FFTGS,
  SGS,
  SPDEGS,

  # ---------
  # LEARNING
  # ---------
  PointwiseLearn

end
