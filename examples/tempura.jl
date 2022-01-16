# build estimators following Toshiya's notes
using Symlens

const ζ⁺ = 1
const ζ⁻ = im

const cϕ = 1    # lensing grad
const cω = -im  # lensing curl
const cα = 1    # rotation
const cτ = 1    # source
const cε = 1    # amplitude

const pϕ = 1
const pω = -1
const pϵ = 1
const ps = 1
const pα = -1

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
    println("building $func (lens)")
    kernels_lens[func] = Symlens.build_l12sum_calculator(
        expr, "lens_$func",
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        evaluate=true
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
    println("building $func (rot)")
    kernels_rot[func] = Symlens.build_l12sum_calculator(
        expr, "rot_$func",
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        evaluate=true
    )
end

#############
# amplitude #
#############
Σ⁰ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)
Σ⁺ = 1/(2ℓ+1)*Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)
Σ⁻ = substitute(Σ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary, which is unsupported for now
Σˣ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)A(ℓ₁)B(ℓ₂)

Γ⁰ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁰(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)
Γ⁺ = 1/(2ℓ+1)*Wₑ⁺(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)
Γ⁻ = substitute(Γ⁺, Dict(ℙ=>-ℙ))  # avoid treating imaginary, which is unsupported for now
Γˣ = 1/(2ℓ+1)*Wₑ⁰(ℓ₁,ℓ,ℓ₂,cε)Wₑ⁺(ℓ₂,ℓ,ℓ₁,cε)A(ℓ₁)B(ℓ₂)

kernels_amp = Dict()
for (expr, func) in zip(
    [Σ⁰, Σ⁺, Σ⁻, Σˣ, Γ⁰, Γ⁺, Γ⁻, Γˣ],
    ["Σ⁰", "Σ⁺", "Σ⁻", "Σˣ", "Γ⁰", "Γ⁺", "Γ⁻", "Γˣ"])
    println("building $func (amp)")
    kernels_amp[func] = Symlens.build_l12sum_calculator(
        expr, "amp_$func",
        Dict(A(ℓ)=>Al, B(ℓ)=>Bl),
        [Al, Bl],
        evaluate=true
    )
end

#-------------------------------------
# Quadratic estimators normalization
#-------------------------------------

function get_kernels_p(est)
    if est == "lens"
        return kernels_lens, pϕ
    elseif est == "rot"
        return kernels_rot, pα
    elseif est == "amp"
        return kernels_amp, pε
    else
        throw
    end
end

function qtt(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al  = @. 1/ocl["TT"]
    Bl  = @. cl["TT"]^2/ocl["TT"]
    res = kernels["Σ⁰"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["TT"]/ocl["TT"]
    res .+= p*kernels["Γ⁰"](lmax, rlmin, rlmax, Al, Al)
    1 ./ res
end

function qte(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al  = @. 1/ocl["TT"]
    Bl  = @. cl["TE"]^2/ocl["EE"]
    res = kernels["Σ⁰"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["TE"]/ocl["TT"]
    Bl .= @. cl["TE"]/ocl["EE"]
    res .+= 2*kernels["Γˣ"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. 1/ocl["EE"]
    Bl .= @. cl["TE"]^2/ocl["TT"]
    res .+= kernels["Σ⁺"](lmax, rlmin, rlmax, Al, Bl)
    1 ./ res
end

function qtb(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al  = @. 1/ocl["BB"]
    Bl  = @. cl["TE"]^2/ocl["TT"]
    res = kernels["Σ⁻"](lmax, rlmin, rlmax, Al, Bl)
    1 ./ res
end

function qee(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al  = @. 1/ocl["EE"]
    Bl  = @. cl["EE"]^2/ocl["EE"]
    res = kernels["Σ⁺"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["EE"]/ocl["EE"]
    res .+= p*kernels["Γ⁺"](lmax, rlmin, rlmax, Al, Al)
    1 ./ res
end

function qbb(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al  = @. 1/ocl["BB"]
    Bl  = @. cl["BB"]^2/ocl["BB"]
    res = kernels["Σ⁺"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["BB"]/ocl["BB"]
    res .+= p*kernels["Γ⁺"](lmax, rlmin, rlmax, Al, Al)
    1 ./ res
end

function qeb(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al = @. 1/ocl["EE"]
    Bl = @. cl["BB"]^2/ocl["BB"]
    res = kernels["Σ⁻"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["BB"]/ocl["BB"]
    Bl .= @. cl["EE"]/ocl["EE"]
    res .+= 2*p*kernels["Γ⁻"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. 1/ocl["BB"]
    Bl .= @. cl["EE"]^2/ocl["EE"]
    res .+= kernels["Σ⁻"](lmax, rlmin, rlmax, Al, Bl)
    1 ./ res
end

function qttte(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al = @. 1/ocl["TT"]
    Bl = @. cl["TT"]*cl["TE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res = kernels["Σ⁰"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["TE"]/ocl["TT"]
    Bl .= @. cl["TT"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res .+= p*kernels["Γˣ"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["TE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl .= @. cl["TT"]/ocl["TT"]
    res .+= p*kernels["Γ⁰"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl .= @. cl["TE"]*cl["TT"]/ocl["TT"]
    res .+= kernels["Σˣ"](lmax, rlmin, rlmax, Al, Bl)
    res
end

function qttee(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al = @. ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl = @. cl["TT"]*cl["EE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res = kernels["Σˣ"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. ocl["TE"]*cl["EE"]/(ocl["TT"]*ocl["EE"])
    Bl .= @. cl["TT"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res .+= p*kernels["Γˣ"](lmax, rlmin, rlmax, Al, Bl)
    res
end

function qteee(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al = @. ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl = @. cl["TE"]*cl["EE"]/ocl["EE"]
    res = kernels["Σˣ"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. cl["TE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl .= @. cl["EE"]/ocl["EE"]
    res .+= p*kernels["Γ⁺"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. (ocl["TE"]*cl["EE"])/(ocl["TT"]*ocl["EE"])
    Bl .= @. (cl["TE"]*ocl["EE"])/(ocl["EE"]^2)
    res .+= p*kernels["Γˣ"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. 1/ocl["EE"]
    Bl .= @. cl["TE"]*cl["EE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res .+= kernels["Σ⁺"](lmax, rlmin, rlmax, Al, Bl)
    res
end

function qtbeb(est, lmax, rlmin, rlmax, cl, ocl)
    kernels, p = get_kernels_p(est)
    Al = @. cl["TE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    Bl = @. cl["BB"]/ocl["BB"]
    res = p*kernels["Γ⁻"](lmax, rlmin, rlmax, Al, Bl)
    Al .= @. 1/ocl["BB"]
    Bl .= @. cl["TE"]*cl["EE"]*ocl["TE"]/(ocl["TT"]*ocl["EE"])
    res .+= kernels["Σ⁻"](lmax, rlmin, rlmax, Al, Bl)
    res
end
