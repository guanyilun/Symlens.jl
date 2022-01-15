module Lensing

using wignerd

function qtt(lmax_p, cltt, nltt)
    lmax_t = length(cltt) - 1

    # initialize gl quadrature
    glq = (lmax_t*2+lmax_p+1)/2 |> round |> Int |> glquad

    # common factors
    ℓ = collect(0:lmax_t)
    llp1      = @. ℓ*(ℓ+1)
    div_dl    = @. 1/(cltt+nltt)
    cl_div_dl = @. cltt/(cltt+nltt)

    # get zeta terms: convert to angular space
    zeta_00   = cf_from_cl(glq, 0,  0, div_dl, prefactor=true)
    zeta_01_p = cf_from_cl(glq, 0,  1, sqrt.(llp1) .* cl_div_dl, prefactor=true)
    zeta_01_m = cf_from_cl(glq, 0, -1, sqrt.(llp1) .* cl_div_dl, prefactor=true)
    zeta_11_p = cf_from_cl(glq, 1,  1, llp1 .* cltt .* cl_div_dl, prefactor=true)
    zeta_11_m = cf_from_cl(glq, 1, -1, llp1 .* cltt .* cl_div_dl, prefactor=true)

    # back to ell space
    nlpp_term_1 = cl_from_cf(glq, -1, -1, lmax_p, zeta_00.*zeta_11_p .- zeta_01_p.^2)
    nlpp_term_2 = cl_from_cf(glq,  1, -1, lmax_p, zeta_00.*zeta_11_m .- zeta_01_p.*zeta_01_m)

    @. (π * llp1 * (nlpp_term_1 + nlpp_term_2)).^-1
end

function qeb(lmax_p, clbb, clee, nleb)
    lmax_e  = length(clee) - 1
    lmax_b  = length(clbb) - 1

    # initialize gl quadrature
    glq = (lmax_t*2+lmax_p+1)/2 |> round |> Int |> glquad

    # common factors
    ℓ         = collect(0:lmax_e)
    cl_en     = @. clee^2 / (clee+nleb)
    cl_en33   = @. cl_en * (ℓ-2) * (ℓ+3)
    cl_en31   = @. cl_en * ((ℓ-1) * (ℓ+2) * (ℓ-2) * (ℓ+3))^0.5
    cl_en11   = @. cl_en * (ℓ-1) * (ℓ+2)

    ℓ = collect(0:lmax_b)
    cl_bn22 = @. (clbb+nleb)^-1

    # convert to angular space
    zeta_en33_p = cf_from_cl(glq, 3,  3, cl_en33, prefactor=true)
    zeta_en33_m = cf_from_cl(glq, 3, -3, cl_en33, prefactor=true)
    zeta_en31_p = cf_from_cl(glq, 3,  1, cl_en31, prefactor=true)
    zeta_en31_m = cf_from_cl(glq, 3, -1, cl_en31, prefactor=true)
    zeta_en11_p = cf_from_cl(glq, 1,  1, cl_en11, prefactor=true)
    zeta_en11_m = cf_from_cl(glq, 1, -1, cl_en11, prefactor=true)

    zeta_bn22_p = cf_from_cl(glq, 2,  2, cl_bn22, prefactor=true)
    zeta_bn22_m = cf_from_cl(glq, 2, -2, cl_bn22, prefactor=true)

    # back to ell space
    nlpp_out_p  = cl_from_cf(glq, 1,  1, lmax_p, @. zeta_en33_p*zeta_bn22_p - 2*zeta_en31_m*zeta_bn22_m + zeta_en11_p*zeta_bn22_p)
    nlpp_out_m  = cl_from_cf(glq, 1, -1, lmax_p, @. zeta_en33_m*zeta_bn22_m - 2*zeta_en31_p*zeta_bn22_p + zeta_en11_m*zeta_bn22_m)

    @. (π/4 * (ℓ*(ℓ+1)) * (nlpp_out_p-nlpp_out_m))^-1
end

end # module
