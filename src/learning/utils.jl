# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LearnedModel(𝓂, θ)

An object that stores a learning model `m`
along with its learned parameters `θ`.
"""
struct LearnedModel{ℳ,Θ}
  𝓂::ℳ
  θ::Θ
end

"""
    learn(𝒯, 𝒟, 𝓂)

Learn the task `𝒯` with geospatial data `𝒟`
using a learning model `𝓂`.
"""
function learn(𝒯::LearningTask, 𝒟, 𝓂)
  # retrieve table of values
  table = values(𝒟)

  # learn model with table
  if issupervised(𝒯)
    X = table |> Select(features(𝒯))
    y = Tables.getcolumn(table, label(𝒯))
    θ, _, __ = MI.fit(𝓂, 0, X, y)
  else
    X = table |> Select(features(𝒯))
    θ, _, __ = MI.fit(𝓂, 0, X)
  end

  # return learned model
  LearnedModel(𝓂, θ)
end

"""
    perform(𝒯, 𝒟, 𝓂̂)

Perform the task `𝒯` with geospatial data `𝒟` using
a *learned* model `𝓂̂`.
"""
function perform(𝒯::LearningTask, 𝒟, 𝓂̂)
  # unpack model and learned parameters
  𝓂, θ = 𝓂̂.𝓂, 𝓂̂.θ

  # retrieve table of values
  table = values(𝒟)

  # apply model to the data
  X = table |> Select(features(𝒯))
  ŷ = MI.predict(𝓂, θ, X)

  # post-process result
  var = outputvars(𝒯)[1]
  val = if issupervised(𝒯)
    isprobabilistic(𝓂) ? mode.(ŷ) : ŷ
  else
    ŷ
  end

  ctor = constructor(typeof(𝒟))
  dom  = domain(𝒟)
  tab  = (; var=>val)
  dat  = Dict(paramdim(dom) => tab)

  ctor(dom, dat)
end
