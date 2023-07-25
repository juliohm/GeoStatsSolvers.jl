elunit(x) = typeunit(eltype(x))
typeunit(::Type) = NoUnits
typeunit(::Type{Q}) where {Q<:Quantity} = unit(Q)

uadjust(::Units, x) = x
function uadjust(U::AffineUnits, x)
  A = absoluteunit(U)
  uconvert.(A, x)
end
