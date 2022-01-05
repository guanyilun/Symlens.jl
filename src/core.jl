using Symbolics, SymbolicUtils, SymbolicUtils.Code

@syms w3j(l1::Number,l2::Number,l3::Number,m1::Number,m2::Number,m3::Number)::Real
@syms F(l1::Number,l2::Number,l3::Number,s::Number)::Real
@syms wigd(l::Number,s1::Number,s2::Number)::Real

rules_F_to_w3j = Vector(undef,2)
rules_F_to_w3j[1] = @rule F(~l1,~l2,~l3,~s) => (-~l1*(~l1+1) + ~l2*(~l2+1) + ~l3*(~l3+1)) * ((2*~l1+1)*(2*~l2+1)*(2*~l3+1))^(1/2) * w3j(~l1,~l2,~l3,-~s,~s,0)  # definition
rules_F_to_w3j[2] = @rule F(~l1,~l2,~l3,~s) => -((2*~l1+1)*(2*~l2+1)*~l3*(~l3+1)*(2*~l3+1))^(1/2) *
    (((~l2-~s)*(~l2+s+1))^(1/2)*w3j(~l1,~l2,~l3,-~s,~s+1,-1) + ((~l2+~s)*(~l2-~s+1))^(1/2)*w3j(~l1,~l2,~l3,-~s,~s-1,1))  # recursive rule

rules_w3j_to_wigd = Vector(undef, 7)
# 3 permutation rules
rules_w3j_to_wigd[1] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l1,~l2,~l3,~s1p,~s2p,~s3p) => wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
rules_w3j_to_wigd[2] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l2,~l3,~l1,~s2p,~s3p,~s1p) => wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
rules_w3j_to_wigd[3] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l3,~l1,~l2,~s3p,~s1p,~s2p) => wigd(~l1,~s1,~s1p)*wigd(~l2,~s2,~s2p)*wigd(~l3,~s3,~s3p)
# 3 flipping rules
rules_w3j_to_wigd[4] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l2,~l1,~l3,~s2p,~s1p,~s3p) => wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
rules_w3j_to_wigd[5] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l3,~l2,~l1,~s3p,~s2p,~s1p) => wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
rules_w3j_to_wigd[6] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)*w3j(~l1,~l3,~l2,~s1p,~s3p,~s2p) => wigd(~l1,~s1,-~s1p)*wigd(~l2,~s2,-~s2p)*wigd(~l3,~s3,-~s3p)
# special case of squaring
rules_w3j_to_wigd[7] = @acrule w3j(~l1,~l2,~l3,~s1,~s2,~s3)^2 => wigd(~l1,~s1,~s1)*wigd(~l2,~s2,~s2)*wigd(~l3,~s3,~s3)

rules_wigd = Vector(undef, 5)
rules_wigd[1] = @rule wigd(~l,~s1::(x->x<0),~s2) => (-1)^(~s1-~s2)*wigd(~l,-~s1,-~s2)
rules_wigd[2] = @rule wigd(~l,~s1,~s2) => ~s1 < ~s2 ? (-1)^(~s1-~s2)*wigd(~l,~s2,~s1) : nothing
rules_wigd[3] = @rule wigd(~l,~s1,~s2) => ~s1 < -~s2 ? wigd(~l,-~s2,-~s1) : nothing

function get_variables_deep(expr)
    vars = Symbolics.get_variables(expr)
    res = []
    for var in vars
        if istree(var)
            for arg in arguments(var)
                push!(res, get_variables_deep(arg)...)
            end
        else
            push!(res, var)
        end
    end
    collect(Set(res))  # drop duplicates
end

function is_atomic(expr)
    vars = get_variables_deep(expr)
    length(vars)>1 ? false : true
end

function are_factors_atomic(expr)
    all(map(ex->is_atomic(ex), arguments(expr)))
end

is_term_of(expr, syms) = (vars = get_variables_deep(expr); Set(vars) == Set(syms))

function factorize_ells(expr)
    # factorize expression to different symbols
    # factors need to be atomic in ells
    @assert are_factors_atomic(expr) "expr is not factorizable"
    op = operation(expr)
    args = arguments(expr)

    syms = get_variables_deep(expr)
    factors = Dict()
    for sym in syms
        factors[sym] = simplify(op(filter(x->is_term_of(x, [sym]), args)...))
    end
    factors
end

count_terms_of(exprs::Vector, target::SymbolicUtils.Sym{SymbolicUtils.FnType{T,N},Nothing}) where {T,N} = (x-> x isa SymbolicUtils.Term ? x.f==target : false).(exprs) |> sum

function get_wigd_s12(expr)
    # get s1 and s2 as a Tuple
    # assume expr is multiply of vectors
    for x ∈ Symbolics.get_variables(expr)
        if x isa SymbolicUtils.Term && x.f == wigd
            return Tuple(arguments(x)[2:end])
        end
    end
    return nothing
end

function drop_wigd(expr)
    # drop wigner d function in the expression, we assume the expression
    # is a product of vector elements or a division whose numerator is a
    # multiplication product.
    @assert (expr isa SymbolicUtils.Mul) || (expr isa SymbolicUtils.Div)
    target = expr isa SymbolicUtils.Mul ? expr : expr.num
    args = filter(x -> !(x isa SymbolicUtils.Term && x.f == wigd), arguments(target))
    new_target = operation(target)(args...)
    expr isa SymbolicUtils.Mul ? new_target : new_target / expr.den
end

# make sure expression of interests is a division
# res = let expr = Nlϕϕ_TT
function build_cl_cf_tables(expr)
    @syms ℓ,ℓ₁,ℓ₂
    @assert expr isa SymbolicUtils.Div
    # make sure denominator is atomic in ells
    @assert are_factors_atomic(expr.den)
    # I. treat denominator, which is easy as it is assumed to be factorizable
    den_factors = factorize_ells(expr.den)
    # II. treat numerator, it will contain wigner 3j related quantities so is more involved
    #  1. if numerator contains F function, use recursive relation to expand it into
    #     wigner 3j symbols
    step1 = simplify(expr.num, RuleSet([rules_F_to_w3j[2]]))
    #  2. expand powers of polynomials if there is any
    step2 = expand(step1)
    #  3. convert wigner 3j products to wigner d matrices
    step3 = simplify(step2, RuleSet(rules_w3j_to_wigd))
    #     make sure no wigner 3j is left
    @assert count_terms_of(Symbolics.get_variables(step3), w3j) == 0 "some wigner 3j left, fail!"
    #  4. massage wigner d to use index convention that the
    #     first `s` is always the larger (in mag) and positive.
    #     Ideally, it will reduce the number of terms
    step4 = simplify(step3, RuleSet(rules_wigd))
    #  5. now if success, we should get an Addition of various factors
    #     each of which is a multiplication with atomic factors.
    num_terms = arguments(step4)
    @assert step4 isa SymbolicUtils.Add "step4 is not of Add type, fail!"
    @assert all([x isa SymbolicUtils.Mul for x ∈ num_terms]) "some factors are not multiplicative products, fail!"
    @assert all([are_factors_atomic(x) for x ∈ num_terms]) "some factors not atomic, fail!"
    #     and make sure we only have three wigner d matrices in each term
    @assert all([count_terms_of(Symbolics.get_variables(term), wigd) == 3 for term ∈ num_terms]) "more than 3 wigd found in some term, fail!"
    #  6. now we are probably safe to proceed to factorize, factorize each term into l, l1, l2 parts
    num_factors_table = factorize_ells.(num_terms)
    #  7. now it's time to incorporate denominator, match l, l1 and l2, modified in place in
    #     num_factors_table
    for num_factors ∈ num_factors_table, s ∈ vars
        num_factors[s] = SymbolicUtils.simplify_div(num_factors[s] / den_factors[s])  # drop common factors
    end
    #  8. collect unique terms for l1 and l2 respectively
    unique_terms_l12 = Dict(sym => collect(Set(map(x->x[sym], num_factors_table))) for sym in [ℓ₁,ℓ₂])  # loop over l1, l2 and collect unique factors
    #  9. since we summe over l1 and l2 eventually and we have factored them out, they can be used interchangably
    #     so we should count them together, this allows us to minimize calculation. To do that we temperarily
    #     change l1 l2 to l
    unique_terms = (x->substitute(x, Dict(ℓ₂ => ℓ₁))).([unique_terms_l12[ℓ₁]...,unique_terms_l12[ℓ₂]...]) |> Set |> collect
    # 10. assign a variable to each of the terms and substitute them back to equation
    @variables ζ[1:length(unique_terms)]
    l1_map = Dict(unique_terms[i] => ζ[i] for i=1:length(unique_terms))
    l2_map = Dict(substitute(unique_terms[i], Dict(ℓ₁=>ℓ₂)) => ζ[i] for i=1:length(unique_terms))
    # 11. now we remerge l,l1,l2 expressions in each term, while keeping track of which (s1, s2)
    #     index it is from wrt wigd matrix.
    cf2cl_table = Dict()
    for num_factors ∈ num_factors_table
        term = reduce(*, [substitute(num_factors[s], Dict(l1_map...,l2_map...)) for s ∈ [ℓ₁,ℓ₂]])
        # get wigner d spins (s1, s2)
        s12  = get_wigd_s12(num_factors[ℓ])
        # get factors apart from wigd
        term *= drop_wigd(num_factors[ℓ])
        # collect all terms with a given (s1, s2) pair
        if haskey(cf2cl_table, s12); cf2cl_table[s12] += term;
        else cf2cl_table[s12] = term end
    end
    cl2cf_table = Dict(v => (get_wigd_s12(k),substitute(drop_wigd(k), Dict(ℓ₂=>ℓ, ℓ₁=>ℓ))) for (k,v) in l1_map)
    # with both of these tables, we should be good to build
    # our function manually
    (cf2cl_table, cl2cf_table)
end

function build_wigd_calls(cl2cf_table, cf2cl_table, rename_table)
    nexprs = length(cl2cf_table)
    exprs = []
    # create n local variables to represent each zeta
    rename = (x) -> substitute(x,
        Dict(k => SymbolicUtils.Sym{Number}(Symbol("zeta_$i"))
            for (i, (k, _)) ∈ enumerate(cl2cf_table)))
    # build assignment expression
    append!(exprs, [:($(rename(k).name) = @__dot__ $(toexpr(substitute(v[2], rename_table))))
                    for (k,v) in cl2cf_table])

    # build cl_from_cl expression, reuse variable names
    append!(exprs, [:($(rename(k).name) = cf_from_cl(glq, $(v[1]...), $(toexpr(rename(k)))))
                    for (k,v) in cl2cf_table])
    # build cf_from_cl expression and add inplace
    for (i,(k,v)) in enumerate(cf2cl_table)
        if i == 1
            push!(exprs, :(res   = cl_from_cf(glq,$(k...),lmax, @__dot__ $(toexpr(rename(v))))))
        else
            push!(exprs, :(res .+= cl_from_cf(glq,$(k...),lmax, @__dot__ $(toexpr(rename(v))))))
        end
    end
    exprs
end

function build_l12sum_calculator(name, cl2cf_table, cf2cl_table, rename_table, args; evaluate=false)
    name = name isa String ? Symbol(name) : name
    f = :(function $(name)(glq, lmax, $(map(x->getfield(x,:name),args)...))
            ℓ = collect(0:max(length.([$(args...)])))
            $(build_wigd_calls(cl2cf_table, cf2cl_table, rename_table)...)
            res
     end)
    evaluate ? eval(f) : f
end