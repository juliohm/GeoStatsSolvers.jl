# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LearnedModel(model, theta)

An object that stores a learning `model`
along with its learned parameters `theta`.
"""
struct LearnedModel{ℳ,Θ}
  𝓂::ℳ
  θ::Θ
end

"""
    learn(task, data, model)

Learn the `task` with geospatial `data`
using a learning `model` from MLJ.jl and
returns a *learned* model.
"""
function learn(𝒯::LearningTask, 𝒟, 𝓂)
  # retrieve table of values
  table = values(𝒟)
  cols = Tables.columns(table)

  # learn model with table
  X = table |> Select(features(𝒯))
  y = Tables.getcolumn(cols, label(𝒯))
  R = MI.reformat(𝓂, X, y)
  θ, _, __ = MI.fit(𝓂, 0, R...)

  # return learned model
  LearnedModel(𝓂, θ)
end

"""
    perform(task, data, lmodel)

Perform the `task` with geospatial `data`
using a *learned* model `lmodel` obtained
from a previous [`learn`](@ref) call.
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
