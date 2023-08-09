# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

elunit(x) = typeunit(eltype(x))

typeunit(::Type) = NoUnits
typeunit(::Type{Q}) where {Q<:Quantity} = unit(Q)

uadjust(x) = uadjust(elunit(x), x)
uadjust(::Units, x) = x
function uadjust(U::AffineUnits, x)
  A = absoluteunit(U)
  uconvert.(A, x)
end

_searchdists!(neighbors, distances, center, metric, domain, method) =
  searchdists!(neighbors, distances, center, method)

function _searchdists!(neighbors, distances, center, metric, domain, method::GlobalSearch)
  nneigh = search!(neighbors, center, method)
  inds = view(neighbors, 1:nneigh)
  xₒ = coordinates(center)
  for i in inds
    x = coordinates(centroid(domain, i))
    distances[i] = metric(xₒ, x)
  end
  nneigh
end
