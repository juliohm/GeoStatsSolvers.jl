# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointwiseLearn(model)

A learning solver that converts geospatial data to a tabular format
with features (and possibly labels) for each point, and then solves
the problem with classical statistical learning `model`.

## Parameters

* `model` - Learning model implementing the `MLJModelInterface.jl`.

## References

* Hoffimann et al. 2020. [Geostatistical Learning: Challenges and Opportunities]
  (https://arxiv.org/abs/2102.08791)
"""
struct PointwiseLearn{M} <: LearningSolver
  model::M
end

function solve(problem::LearningProblem, solver::PointwiseLearn)
  sdata = sourcedata(problem)
  tdata = targetdata(problem)
  ptask = task(problem)
  model = solver.model

  # assert model is compatible with task
  @assert iscompatible(model, ptask) "$model is not compatible with $ptask"

  # learn model on source data
  lmodel = learn(ptask, sdata, model)

  # apply model to target data
  perform(ptask, tdata, lmodel)
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::PointwiseLearn)
  print(io, "PointwiseLearn")
end

function Base.show(io::IO, ::MIME"text/plain", solver::PointwiseLearn)
  println(io, solver)
  print(io, "  └─model ⇨ ")
  show(IOContext(io, :compact => true), solver.model)
end
