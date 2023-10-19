# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    PointwiseLearn(model)

The pointwise learning solver introduced by Hoffimann 2021 to
illustrate the challenges of applying classical machine learning
"point-by-point" with geospatial data.

## Parameters

* `model` - Learning model supported by the StatsLearnModels.jl package.

### References

* Hoffimann et al. 2021. [Geostatistical Learning: Challenges and Opportunities]
  (https://www.frontiersin.org/articles/10.3389/fams.2021.689393/full)

### Notes

* This solver is equivalent to a vanilla framework for
  machine learning (e.g. MLJ.jl). It exists for didactic
  and practical purposes since it automates the process
  of learning functions over geospatial domains.
"""
struct PointwiseLearn{M} <: LearningSolver
  model::M
end

function solve(problem::LearningProblem, solver::PointwiseLearn)
  sdata = sourcedata(problem)
  tdata = targetdata(problem)
  ptask = task(problem)
  model = solver.model

  # learn model on source data
  lmodel = Learn(values(sdata), model, features(ptask) => label(ptask))

  # apply model to target data
  pred = lmodel(values(tdata))
  georef(pred, domain(tdata))
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
