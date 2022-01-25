module Symlens

using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using SymbolicUtils.Rewriters
using SyntaxTree  # for better readability of returned expression
using Wignerd

# convenience
export @syms, substitute

# reserved symbolic names
# Fₒₖ represents F function in Okamoto & Hu lensing paper
# ̱̱ℙ represents parity factor, i.e., (-1)^(ℓ+ℓ₁+ℓ₂)
@syms w3j(l1::Number,l2::Number,l3::Number,m1::Number,m2::Number,m3::Number)::Real
@syms Fₒₖ(l1::Number,l2::Number,l3::Number,s::Number)::Real
@syms wigd(l::Number,s1::Number,s2::Number)::Real
@syms ℚ(l::Number) ℙ::Number  # \bbQ \bbP, to avoid conflicts with common variable names
@syms ℓ ℓ₁ ℓ₂

export w3j, Fₒₖ, wigd, ℓ, ℓ₁, ℓ₂, ℙ

include("utils.jl")
include("rewriters.jl")
include("core.jl")
export build_l12sum_calculator, build_cl_cf_tables, factorize_wigd

include("Prebuilt/Prebuilt.jl")

end
