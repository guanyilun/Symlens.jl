# Fₒₖ is the ₛF function given in Okamoto & Hu paper, hence the name
rules_F_to_w3j = Vector(undef,2)
# definition
rules_F_to_w3j[1] = @rule Fₒₖ(~l1,~l2,~l3,~s) => (-~l1*(~l1+1) + ~l2*(~l2+1) + ~l3*(~l3+1)) * ((2*~l1+1)*(2*~l2+1)*(2*~l3+1)/(16π))^(1/2) * w3j(~l1,~l2,~l3,-~s,~s,0)
# recursion
rules_F_to_w3j[2] = @rule Fₒₖ(~l1,~l2,~l3,~s) => -((2*~l1+1)*(2*~l2+1)*~l3*(~l3+1)*(2*~l3+1)/(16π))^(1/2) * (((~l2-~s)*(~l2+s+1))^(1/2)*w3j(~l1,~l2,~l3,-~s,~s+1,-1) + ((~l2+~s)*(~l2-~s+1))^(1/2)*w3j(~l1,~l2,~l3,-~s,~s-1,1))

rules_w3j_to_wigd = Vector(undef, 7)
# 3 permutation rules
rules_w3j_to_wigd[1] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l1,~l2,~l3,~s1p,~s2p,~s3p) => 0.5*wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
rules_w3j_to_wigd[2] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l2,~l3,~l1,~s2p,~s3p,~s1p) => 0.5*wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
rules_w3j_to_wigd[3] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l3,~l1,~l2,~s3p,~s1p,~s2p) => 0.5*wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
# 3 flipping rules
rules_w3j_to_wigd[4] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l2,~l1,~l3,~s2p,~s1p,~s3p) => 0.5*wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
rules_w3j_to_wigd[5] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l3,~l2,~l1,~s3p,~s2p,~s1p) => 0.5*wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
rules_w3j_to_wigd[6] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l1,~l3,~l2,~s1p,~s3p,~s2p) => 0.5*wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
# special case of squaring
rules_w3j_to_wigd[7] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)^2 => 0.5*wigd(~l1,~s1,~s1)*wigd(~l2,~s2,~s2)*wigd(~l3,~s3,~s3)

# s1, s2 convention
rules_wigd = Vector(undef, 3)
rules_wigd[1] = @rule wigd(~l,~s1::(x->x<0),~s2) => (-1)^(~s1-~s2)*wigd(~l,-~s1,-~s2)
rules_wigd[2] = @rule wigd(~l,~s1,~s2) => ~s1 < ~s2 ? (-1)^(~s1-~s2)*wigd(~l,~s2,~s1) : nothing
rules_wigd[3] = @rule wigd(~l,~s1,~s2) => ~s1 < -~s2 ? wigd(~l,-~s2,-~s1) : nothing

# Q(l) = (-1)^l. The rules below assumes that Q comes entirely from
# (-1)^(l1+l2+l3), and integrated over \theta, which is the only case
# the transformation rule below is true
rules_Q = Vector(undef, 2)
rules_Q[1] = @acrule ℚ(~l)^(~n)*wigd(~l,~s1,~s2) => wigd(~l,~s1,(-1)^(~n)*~s2)
rules_Q[2] = @acrule ℚ(~l)*wigd(~l,~s1,~s2) => wigd(~l,~s1,-1*~s2)

P_to_Q(expr) = substitute(expr, Dict(ℙ=>ℚ(ℓ₁)*ℚ(ℓ₂)*ℚ(ℓ)))
F_to_w3j(expr) = rules_F_to_w3j[2](expr)
w3j_to_wigd(expr) = Fixpoint(Postwalk(Chain(rules_w3j_to_wigd)))(expr)
resolve_Q(expr) = Fixpoint(Postwalk(Chain(rules_Q)))(expr)
wigd_convention(expr) = Postwalk(Chain(rules_wigd))(expr)

function expand_nonatomic(expr)
    args = arguments(expr)
    op = operation(expr)
    atomics = filter(is_atomic, args)
    # represent each atomic factor with a temperary variable
    # then expand and change back
    forward = Dict(x => SymbolicUtils.Sym{Number}(gensym()) for x in atomics)
    backward = reverse(forward)
    substitute(expr, forward) |> expand |> x->substitute(x, backward)
end
