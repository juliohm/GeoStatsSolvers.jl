# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoStatsSolvers

using Meshes
using GeoTables
using GeoStatsBase
using Variography
using GeoStatsModels

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

# aliases
import MLJModelInterface
const MI = MLJModelInterface

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
include("learning/utils.jl")

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

  # meta solvers
  CookieCutter,

  # ---------
  # LEARNING
  # ---------
  PointwiseLearn,

  # utilities
  learn,
  perform

end
