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

  if isnothing(maxneighbors) || maxneighbors == nelem
    # global search with all elements
    GlobalSearch(domain)
  else
    if maxneighbors > nelem
      throw(ArgumentError("maxneighbors must be smaller than number of elements"))
    end

    if isnothing(neighborhood)
      # nearest neighbor search with a metric
      KNearestSearch(domain, maxneighbors, metric=metric)
    else
      # neighbor search with ball neighborhood
      KBallSearch(domain, maxneighbors, neighborhood)
    end
  end
end

"""
    kriging_ui(domain, variogram, mean, degree, drifts)

Return the appropriate Kriging estimator for the `domain` based on
end-user inputs such as `variogram`, `mean`, `degree` and `drifts`.
"""
function kriging_ui(domain, variogram, mean, degree, drifts)
  if drifts ≠ nothing
    ExternalDriftKriging(variogram, drifts)
  elseif degree ≠ nothing
    UniversalKriging(variogram, degree, embeddim(domain))
  elseif mean ≠ nothing
    SimpleKriging(variogram, mean)
  else
    OrdinaryKriging(variogram)
  end
end
