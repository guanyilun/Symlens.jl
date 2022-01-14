module symlens

using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using SymbolicUtils.Rewriters

using wignerd

# reserved symbolic names
@syms w3j(l1::Number,l2::Number,l3::Number,m1::Number,m2::Number,m3::Number)::Real
@syms Fₒₖ(l1::Number,l2::Number,l3::Number,s::Number)::Real
@syms wigd(l::Number,s1::Number,s2::Number)::Real
@syms ℚ(l::Number) ℙ::Number  # \bbQ \bbP, to avoid conflicts with common variable names
@syms ℓ ℓ₁ ℓ₂

export w3j, Fₒₖ, wigd, ℓ, ℓ₁, ℓ₂, ℙ

include("rewriters.jl")
include("utils.jl")
include("core.jl")

end
