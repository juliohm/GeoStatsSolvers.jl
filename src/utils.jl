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
