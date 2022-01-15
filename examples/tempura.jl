# build estimators following Toshiya's notes
using Symlens

ζ⁺ = 1
ζ⁻ = im

cϕ = 1    # lensing grad
cω = -im  # lensing curl
cα = 1    # rotation
cτ = 1    # source
cε = 1    # amplitude

γ(ℓ₁,ℓ₂,ℓ₃) = ((2ℓ₁+1)*(2ℓ₂+1)*(2ℓ₃+1)/4π)^(1/2)

# parity operator: ̱ℙ ≡ (-1)^(\sum_i l_i)
q⁺(c) = c*(1+ℙ)/2.0
q⁻(c) = c*(1-ℙ)/2.0
qₓ⁺(c) = c*(1+c^2*ℙ)/2.0
qₓ⁻(c) = c*(1-c^2*ℙ)/2.0

a(ℓ,s) = -((ℓ-s)*(ℓ+s+1)/2.0)^(1/2)
a⁺(ℓ) = a(ℓ,2)
a⁻(ℓ) = a(ℓ,-2)

# lensing
Wₓ⁰(ℓ₁,ℓ₂,ℓ₃,c) = -2*a(ℓ₂,0)*a(ℓ₃,0)*qₓ⁺(c)*γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,0,1,-1)
Wₓ⁺(ℓ₁,ℓ₂,ℓ₃,c) = -ζ⁺*qₓ⁺(c)*γ(ℓ₁,ℓ₂,ℓ₃)*a(ℓ₂,0)*(a⁺(ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,1,-3)+c^2*a⁻(ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,-1,-1))
Wₓ⁻(ℓ₁,ℓ₂,ℓ₃,c) = -ζ⁻*qₓ⁻(c)*γ(ℓ₁,ℓ₂,ℓ₃)*a(ℓ₂,0)*(a⁺(ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,1,-3)+c^2*a⁻(ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,-1,-1))

# rotation
Wₐ⁺(ℓ₁,ℓ₂,ℓ₃,c) = 2im*ζ⁺*q⁻(c)*γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,0,-2)
Wₐ⁻(ℓ₁,ℓ₂,ℓ₃,c) = 2im*ζ⁻*q⁺(c)*γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,0,-2)

# amplitude
Wₑ⁰(ℓ₁,ℓ₂,ℓ₃,c) = γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,0,0,0)
Wₑ⁺(ℓ₁,ℓ₂,ℓ₃,c) = ζ⁺*q⁺(c)*γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,0,-2)
Wₑ⁻(ℓ₁,ℓ₂,ℓ₃,c) = ζ⁻*q⁻(c)*γ(ℓ₁,ℓ₂,ℓ₃)*w3j(ℓ₁,ℓ₂,ℓ₃,2,0,-2)

#-----------------------------------
# normalization kernel calculation 
#-----------------------------------

# general definition:
# each term below is missing A(ℓ₁)B(ℓ₂) and a summation over ℓ₁ and ℓ₂
# Σ⁰(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁰(ℓ₁,ℓ,ℓ₂,c))Wₓ⁰(ℓ₁,ℓ,ℓ₂,c)
# Γ⁰(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁰(ℓ₁,ℓ,ℓ₂,c))Wₓ⁰(ℓ₂,ℓ,ℓ₁,c)

# Σ⁺(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁺(ℓ₁,ℓ,ℓ₂,c))Wₓ⁺(ℓ₁,ℓ,ℓ₂,c)
# Γ⁺(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁺(ℓ₁,ℓ,ℓ₂,c))Wₓ⁺(ℓ₂,ℓ,ℓ₁,c)

# Σ⁻(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁻(ℓ₁,ℓ,ℓ₂,c))Wₓ⁻(ℓ₁,ℓ,ℓ₂,c)
# Γ⁻(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*conj(Wₓ⁻(ℓ₁,ℓ,ℓ₂,c))Wₓ⁻(ℓ₂,ℓ,ℓ₁,c)

# Σˣ(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,c)Wₓ⁺(ℓ₁,ℓ,ℓ₂,c)
# Γˣ(ℓ,ℓ₁,ℓ₂,c) = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,c)Wₓ⁺(ℓ₂,ℓ,ℓ₁,c)

@syms A(ℓ) B(ℓ) ℓ ℓ₁ ℓ₂
@syms Al Bl

###########
# lensing #
###########
Σ⁰ = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁰(ℓ₁,ℓ,ℓ₂,cϕ)A(ℓ₁)B(ℓ₂)
Σ⁺ = 1/(2ℓ+1)*Wₓ⁺(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁺(ℓ₁,ℓ,ℓ₂,cϕ)A(ℓ₁)B(ℓ₂)
Σ⁻ = substitute(Σ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary, which is unsupported for now
Σˣ = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁺(ℓ₁,ℓ,ℓ₂,cϕ)A(ℓ₁)B(ℓ₂)
# l₂ <-> l₁ in second W
Γ⁰ = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁰(ℓ₂,ℓ,ℓ₁,cϕ)A(ℓ₁)B(ℓ₂) 
Γ⁺ = 1/(2ℓ+1)*Wₓ⁺(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁺(ℓ₂,ℓ,ℓ₁,cϕ)A(ℓ₁)B(ℓ₂)
Γ⁻ = substitute(Γ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary
Γˣ = 1/(2ℓ+1)*Wₓ⁰(ℓ₁,ℓ,ℓ₂,cϕ)Wₓ⁺(ℓ₂,ℓ,ℓ₁,cϕ)A(ℓ₁)B(ℓ₂)

kernels_lens = Dict()
for (expr, func) in zip(
    [Σ⁰, Σ⁺, Σ⁻, Σˣ, Γ⁰, Γ⁺, Γ⁻, Γˣ],
    ["Σ⁰", "Σ⁺", "Σ⁻", "Σˣ", "Γ⁰", "Γ⁺", "Γ⁻", "Γˣ"])
    println("building $func")
    kernels_lens[func] = Symlens.build_l12sum_calculator(
        expr, func,
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        prefactor=true, evaluate=true
    )
end

############
# rotation #
############

Σ⁺ = 4/(2ℓ+1)*q⁻(cα)*γ(ℓ₁,ℓ,ℓ₂)*w3j(ℓ₁,ℓ,ℓ₂,2,0,-2)*q⁻(cα)*γ(ℓ₁,ℓ,ℓ₂)*w3j(ℓ₁,ℓ,ℓ₂,2,0,-2)
Σ⁻ = substitute(Σ⁺, Dict(ℙ=>-ℙ))
Γ⁺ = 4/(2ℓ+1)*q⁻(cα)*γ(ℓ₁,ℓ,ℓ₂)*w3j(ℓ₁,ℓ,ℓ₂,2,0,-2)*q⁻(cα)*γ(ℓ₂,ℓ,ℓ₁)*w3j(ℓ₂,ℓ,ℓ₁,2,0,-2)
Γ⁻ = substitute(Γ⁺, Dict(ℙ=>-ℙ))

@syms Al Bl
kernels_rot = Dict()
for (expr, func) in zip(
    [Σ⁺, Σ⁻, Γ⁺, Γ⁻],
    ["Σ⁺", "Σ⁻", "Γ⁺", "Γ⁻"])
    println("building $func")
    kernels_rot[func] = Symlens.build_l12sum_calculator(
        expr, func,
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        prefactor=true, evaluate=true
    )
end

#############
# amplitude #
#############
Σ⁰ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)
Σ⁺ = 1/(2ℓ+1)*Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)
Σ⁻ = substitute(Σ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary, which is unsupported for now
Σˣ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)

Σ⁰ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁰(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)
Σ⁺ = 1/(2ℓ+1)*Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)
Σ⁻ = substitute(Σ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary, which is unsupported for now
Σˣ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)

kernels_amp = Dict()
for (expr, func) in zip(
    [Σ⁰, Σ⁺, Σ⁻, Σˣ, Γ⁰, Γ⁺, Γ⁻, Γˣ],
    ["Σ⁰", "Σ⁺", "Σ⁻", "Σˣ", "Γ⁰", "Γ⁺", "Γ⁻", "Γˣ"])
    println("building $func")
    kernels_amp[func] = Symlens.build_l12sum_calculator(
        expr, func,
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        prefactor=true, evaluate=true
    )
end
