# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    searcher_ui(domain, maxneighbors, metric, neighborhood)

Return the appropriate search method over the `domain` based on
end-user inputs such as `maxneighbors`, `metric` and `neighborhood`.
"""
function searcher_ui(domain, maxneighbors, metric, neighborhood)
  # number of domain elements
  nelem = nelements(domain)

  # number of neighbors
  nmax = if isnothing(maxneighbors)
    nelem
  elseif maxneighbors < 1 || maxneighbors > nelem
    @warn "Invalid maximum number of neighbors. Using all elements..."
    nelem
  else
    maxneighbors
  end

  if isnothing(neighborhood)
    # nearest neighbor search with a metric
    KNearestSearch(domain, nmax; metric)
  else
    # neighbor search with ball neighborhood
    KBallSearch(domain, nmax, neighborhood)
  end
end

"""
    kriging_ui(domain, variogram, mean, degree, drifts)

Return the appropriate Kriging estimator for the `domain` based on
end-user inputs such as `variogram`, `mean`, `degree` and `drifts`.
"""
function kriging_ui(domain, variogram, mean, degree, drifts)
  if !isnothing(drifts)
    ExternalDriftKriging(variogram, drifts)
  elseif !isnothing(degree)
    UniversalKriging(variogram, degree, embeddim(domain))
  elseif !isnothing(mean)
    SimpleKriging(variogram, mean)
  else
    OrdinaryKriging(variogram)
  end
end
