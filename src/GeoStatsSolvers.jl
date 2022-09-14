# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoStatsSolvers

using Meshes
using GeoStatsBase
using Variography
using KrigingEstimators

using Tables
using Distributions
using TableTransforms
using NearestNeighbors
using Distances: Euclidean
using Bessels: gamma
using CpuId
using FFTW

using LinearAlgebra
using Statistics
using Random

# aliases
import MLJModelInterface
const MI = MLJModelInterface

import GeoStatsBase: preprocess, solve, solvesingle

include("estimation/idw.jl")
include("estimation/lwr.jl")
include("estimation/kriging.jl")

include("simulation/lu.jl")
include("simulation/fft.jl")
include("simulation/seq.jl")
include("simulation/sgs.jl")
include("simulation/spde.jl")
include("simulation/cookie.jl")

include("learning/pointwise.jl")
include("learning/utils.jl")

export
  # -----------
  # ESTIMATION
  # -----------
  IDW,
  LWR,
  Kriging,

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

  # meta solvers
  CookieCutter,

  # ---------
  # LEARNING
  # ---------
  PointwiseLearn,

  # utilities
  learn, perform

end
