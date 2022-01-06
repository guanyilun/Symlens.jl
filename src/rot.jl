module rot

using wignerd

"""
    clbb_from_claa(lmax_b, clee, claa)

Calculate ClBB from Faraday's rotation based on an input rotational
power spectrum. It is assumed that the power spectra input starts from ell=0

Parameters
- lmax_b: lmax for ClBB
- clee: ClEE power spectrum
- claa: Cl^\alpha\alpha power spectrum (rotational power spectrum)

Returns
- clbb
"""
function clbb_from_claa(lmax, clee, claa)
    npoints = (3*max([lmax, length.([clee, claa])...]...)+1)/2 |> Int
    glq = wignerd.glquad(npoints)

    ls = collect(0:length(claa)-1)
    zeta_00 = cf_from_cl(glq, 0, 0, claa, prefactor=true)

    ls = collect(0:length(clee)-1)
    zeta_m2m2 = cf_from_cl(glq, -2, -2, clee, prefactor=true)
    zeta_p2m2 = cf_from_cl(glq,  2, -2, clee, prefactor=true)
    clbb = (cl_from_cf(glq,  2, 2, lmax, @. zeta_00 * zeta_m2m2) .+
            cl_from_cf(glq, -2, 2, lmax, @. zeta_00 * zeta_p2m2)) .* (4π)
    clbb
end

end  # module
