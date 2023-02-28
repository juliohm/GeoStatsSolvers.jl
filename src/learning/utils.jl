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
  cols  = Tables.columns(table)

  # learn model with table
  X = table |> Select(features(𝒯))
  y = Tables.getcolumn(cols, label(𝒯))
  R = MI.reformat(𝓂, X, y)
  θ, _, __ = MI.fit(𝓂, 0, R...)

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
  R = MI.reformat(𝓂, X)
  ŷ = if isprobabilistic(𝓂)
    MI.predict_mode(𝓂, θ, R...)
  else
    MI.predict(𝓂, θ, R...)
  end

  georef((; label(𝒯) => ŷ), domain(𝒟))
end
