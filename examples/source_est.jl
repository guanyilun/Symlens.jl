using Symlens

@syms y(ℓ)

W = ((2ℓ+1)*(2ℓ₁+1)*(2ℓ₂+1)/4π)^(1/2)*w3j(ℓ₁,ℓ₂,ℓ,0,0,0)*y(ℓ₁)*y(ℓ₂)

@syms f(ℓ)

R = (W^2*f(ℓ₁)*f(ℓ₂))/(2*(2*ℓ+1))

@syms fl yl
Symlens.build_l12sum_calculator(R, "source_est",
                                Dict(f(ℓ)=>fl, y(ℓ)=>yl),
                                [fl, yl]; prefactor=true,
                                post = [:(res .^= -1)],
                                evaluate=false) |> println
