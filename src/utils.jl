elunit(x) = typeunit(eltype(x))
typeunit(::Type) = NoUnits
typeunit(::Type{Q}) where {Q<:Quantity} = unit(Q)

uadjust(::Unitful.Units, x) = x
function uadjust(U::Unitful.AffineUnits, x)
  A = absoluteunit(U)
  uconvert.(A, x)
end

urevert(::Unitful.Units, x) = x
urevert(U::Unitful.AffineUnits, x) = uconvert.(U, x)
