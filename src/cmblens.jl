module cmblens

using wignerd

function qtt(lmax_p, cltt, nltt)
    lmax_t = len(cltt) - 1
    npoints = convert(Integer, (lmax_t*2+lmax_p)/2+1)
    # initialize gl quadrature
    glq = wignerd.glquad(int(lmax_t*2+lmax_p)/2+1)

    # common factors
    ℓ = collect(0:lmax_t)
    prefactor = @. (2*ℓ+1)/(4*π)
    llp1      = @. ℓ*(ℓ+1)
    div_dl    = @. 1/(cltt+nltt)
    cl_div_dl = @. cltt/(cltt+nltt)

    # get zeta terms
    zeta_00   = wignerd.cf_from_cl(glq, 0,  0, prefactor .* div_dl)
    zeta_01_p = wignerd.cf_from_cl(glq, 0,  1, prefactor .* sqrt.(llp1) .* cl_div_dl)
    zeta_01_m = wignerd.cf_from_cl(glq, 0, -1, prefactor .* sqrt.(llp1) .* cl_div_dl)
    zeta_11_p = wignerd.cf_from_cl(glq, 1,  1, prefactor .* llp1 .* cltt .* cl_div_dl)
    zeta_11_m = wignerd.cf_from_cl(glq, 1, -1, prefactor .* llp1 .* cltt .* cl_div_dl)

    nlpp_term_1 = wignerd.cl_from_cf(glq, -1, -1, lmax_p, zeta_00.*zeta_11_p .- zeta_01_p.^2)
    nlpp_term_2 = wignerd.cl_from_cf(glq,  1, -1, lmax_p, zeta_00.*zeta_11_m .- zeta_01_p.*zeta_01_m)

    @. (π * llp1 * (nlpp_term_1 + nlpp_term_2)).^-1
end

function qeb(lmax_p, clbb, clee, nleb)
    lmax_e  = length(clee) - 1
    lmax_b  = length(clbb) - 1
    npoints = convert(Integer, (lmax_p+lmax_b+lmax_e)*0.5 + 1)
    gl = glquad(npoints)

    ℓ         = collect(0:lmax_e)
    prefactor = @. (2*ℓ+1)/(4*π)
    cl_en     = @. clee^2 / (clee+nleb)
    cl_en33   = @. prefactor * cl_en * (ℓ-2) * (ℓ+3)
    cl_en31   = @. prefactor * cl_en * ((ℓ-1) * (ℓ+2) * (ℓ-2) * (ℓ+3))^0.5
    cl_en11   = @. prefactor * cl_en * (ℓ-1) * (ℓ+2)

    ℓ = collect(0:lmax_b)
    cl_bn22 = @. prefactor * (clbb+nleb)^-1

    zeta_en33_p = wignerd.cf_from_cl(glq, 3,  3, cl_en33)
    zeta_en33_m = wignerd.cf_from_cl(glq, 3, -3, cl_en33)
    zeta_en31_p = wignerd.cf_from_cl(glq, 3,  1, cl_en31)
    zeta_en31_m = wignerd.cf_from_cl(glq, 3, -1, cl_en31)
    zeta_en11_p = wignerd.cf_from_cl(glq, 1,  1, cl_en11)
    zeta_en11_m = wignerd.cf_from_cl(glq, 1, -1, cl_en11)

    zeta_bn22_p = wignerd.cf_from_cl(glq, 2,  2, cl_bn22)
    zeta_bn22_m = wignerd.cf_from_cl(glq, 2, -2, cl_bn22)

    nlpp_out_p  = wignerd.cl_from_cf(glq, 1,  1, lmax_p, @. zeta_en33_p*zeta_bn22_p - 2*zeta_en31_m*zeta_bn22_m + zeta_en11_p*zeta_bn22_p)
    nlpp_out_m  = wignerd.cl_from_cf(glq, 1, -1, lmax_p, @. zeta_en33_m*zeta_bn22_m - 2*zeta_en31_p*zeta_bn22_p + zeta_en11_m*zeta_bn22_m)

    @. (π/4 * (ℓ*(ℓ+1)) * (nlpp_out_p-nlpp_out_m))^-1
end

end # module
