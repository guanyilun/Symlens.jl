function factorize_ells(expr)
    # factorize expression to different symbols
    # factors need to be atomic in ells
    @assert are_factors_atomic(expr) "expr is not factorizable"
    op   = operation(expr)
    args = arguments(expr)
    syms = [ℓ, ℓ₁, ℓ₂]
    # we assume it is a multiplication of different atomic factors,
    # like (ℓ+1)(ℓ₂+1).... The only exception allowed here is when the
    # expression is an addition which depends only on one variable,
    # such as 2l+1. Below we will handle this special case
    factors = Dict()
    if expr isa SymbolicUtils.Add
        vars = get_variables_deep(expr)
        @assert length(vars) == 1 "Addition of different variables, not sure how to factor, fail!"
        # now we assume that the entire expr is a factor of one variable only
        factors[vars[1]] = expr
        for s in filter(x->(x !== vars[1]), syms)
            factors[s] = 1
        end
        factors[:const] = 1
    else
        @assert expr isa SymbolicUtils.Mul "Not add nor mul, not sure how to factor, fail!"
        for sym in syms
            terms = filter(x->is_term_of(x, [sym]), args)
            factors[sym] = length(terms) > 0 ? simplify(op(terms...)) : 1
        end
        # don't forget about constant
        constants = filter(x->is_term_of(x, Any[]), args)
        if length(constants)>0; factors[:const] = op(constants...)
        else factors[:const] = 1 end
    end
    factors
end

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
    if (expr isa SymbolicUtils.Term && expr.f == wigd); return 1 end
    @assert (expr isa SymbolicUtils.Mul) || (expr isa SymbolicUtils.Div)
    target = expr isa SymbolicUtils.Mul ? expr : expr.num
    # target could be wigd itself
    if (target isa SymbolicUtils.Term && target.f == wigd)
        new_target = 1
    else
        # otherwise assume it's a multiplication
        @assert target isa SymbolicUtils.Mul "unhandled case"
        args = filter(x -> !(x isa SymbolicUtils.Term && x.f == wigd), arguments(target))
        new_target = operation(target)(args...)
    end
    expr isa SymbolicUtils.Div ? new_target / expr.den : new_target
end

"""try to reduce the number of terms in cf_tables
by grouping similar terms"""
function reduce_cl_table(cl_table)
    reduce_table = Dict()
    for v in cl_table
        k = (v[1],v[3])
        if haskey(reduce_table, k); reduce_table[k] += v[2]
        else reduce_table[k] = v[2] end
    end
    # map back
    [[k[1],v,k[2]] for (k,v) in reduce_table]
end

"""build some factorization tables which contain all information needed to build a function.
prefactor assumes the cf_from_cl calls will have prefactor option turned on

"""
function build_cl_cf_tables(expr; prefactor=false, verbose=false)
    @syms ℓ ℓ₁ ℓ₂
    @assert expr isa SymbolicUtils.Div || expr isa SymbolicUtils.Mul
    isdiv = expr isa SymbolicUtils.Div

    # treat prefactor in different cases
    prefactor && isdiv  && (expr = (16π^2*expr.num) / (expr.den*(2ℓ₁+1)*(2ℓ₂+1));)
    prefactor && !isdiv && (expr = (16π^2*expr) / ((2ℓ₁+1)*(2ℓ₂+1)); isdiv=true)

    # 1. treat denominator, which is easy as it is assumed to be
    # factorizable We first make sure denominator is atomic in ells,
    # and then factorize it into l, l1, l2, const terms.
    if isdiv
        @assert are_factors_atomic(expr.den)
        den_factors = factorize_ells(expr.den)
        verbose && (@show den_factors)
    end
    # 2. treat numerator, it will contain wigner 3j related quantities
    # so is more involved. If we didn't get a division, treat it as
    # the numerator)
    num = isdiv ? expr.num : expr
    step = factorize_wigd(num)
    verbose && (@show step)
    # 3. now if success, we should get an Addition of various factors
    # each of which is a multiplication with atomic factors.
    num_terms = step isa SymbolicUtils.Add ? argument(step) : [step]
    @assert all([x isa SymbolicUtils.Mul for x ∈ num_terms]) "some factors are not multiplicative products, fail!"
    @assert all([are_factors_atomic(x) for x ∈ num_terms]) "some factors not atomic, fail!"
    # and make sure we only have three wigner d matrices in each term
    @assert all([count_terms_of(Symbolics.get_variables(term), wigd) == 3
                 for term ∈ num_terms]) "more than 3 wigd found in some term, fail!"
    # 4. now we are probably safe to proceed to factorize, factorize
    # each term into l, l1, l2 parts
    num_factors_table = factorize_ells.(num_terms)
    verbose && (@show num_factors_table)
    # 5. now it's time to incorporate denominator, if any, to match l,
    # l1 and l2, modified in place in num_factors_table
    if isdiv
        for num_factors ∈ num_factors_table
            for s ∈ [ℓ, ℓ₁, ℓ₂]
                num_factors[s] = simplify_fractions(num_factors[s] / den_factors[s])
            end
            num_factors[:const] = num_factors[:const] / den_factors[:const]  # no need to simplify constants
        end
    end
    verbose && (println("num_factors_table (after div): ", num_factors_table))
    # 6. collect unique terms for l1 and l2 respectively
    unique_terms_l12 = Dict(sym => collect(Set(map(x->x[sym], num_factors_table))) for sym in [ℓ₁,ℓ₂])  # loop over l1, l2 and collect unique factors
    verbose && (@show unique_terms_l12)
    # 7. since we summe over l1 and l2 eventually and we have factored
    # them out, they can be used interchangably so we should count
    # them together, this allows us to minimize calculation. To do
    # that we temperarily change l1 l2 to l
    unique_terms = (x->substitute(x, Dict(ℓ₂ => ℓ₁))).([unique_terms_l12[ℓ₁]...,unique_terms_l12[ℓ₂]...]) |> Set |> collect
    verbose && (@show unique_terms)
    # 8. assign a variable to each of the terms and substitute them back to equation
    @variables ζ[1:length(unique_terms)]
    l1_map = Dict(unique_terms[i] => ζ[i] for i=1:length(unique_terms))
    verbose && (@show l1_map)
    l2_map = Dict(substitute(unique_terms[i], Dict(ℓ₁=>ℓ₂)) => ζ[i] for i=1:length(unique_terms))
    verbose && (@show l2_map)
    # 9. now we remerge l,l1,l2 expressions in each term, while
    # keeping track of which (s1, s2) index it is from wrt wigd
    # matrix.
    cl_table = []
    for num_factors ∈ num_factors_table
        term = reduce(*, [substitute(num_factors[s], Dict(l1_map...,l2_map...)) for s ∈ [ℓ₁,ℓ₂]]) * num_factors[:const]
        # get wigner d spins (s1, s2)
        s12  = get_wigd_s12(num_factors[ℓ])
        # get l terms apart from wigd
        lterm = drop_wigd(num_factors[ℓ])
        push!(cl_table, [s12, term, lterm])
    end
    # try to reduce the number of terms in cl_table by grouping common
    # factors with both of these tables, we should be good to build
    # our function manually
    verbose && (@show cl_table)
    cl_table = reduce_cl_table(cl_table)
    verbose && (println("cl_table (after reduce): ", cl_table))
    cf_table = Dict(v => (get_wigd_s12(k),substitute(drop_wigd(k), Dict(ℓ₁=>ℓ))) for (k,v) in l1_map)
    verbose && (@show cf_table)

    (cl_table, cf_table)
end

function build_wigd_calls(cl_table, cf_table, rename_table; prefactor=false)
    nexprs = length(cf_table)
    exprs = []
    # create n local variables to represent each zeta note that this
    # is not elegant as often times some of these zeta variables are
    # identical, so I should add another unique and variable remapping
    # here at some point <- FIXME
    rename = (x) -> substitute(x,
        Dict(k => SymbolicUtils.Sym{Number}(Symbol("zeta_$i"))
            for (i, (k, _)) ∈ enumerate(cf_table)))
    # build assignment expression
    append!(exprs, [:($(rename(k).name) = @__dot__ $(substitute(v[2], rename_table)))
                    for (k,v) in cf_table])
    # build cl_from_cl expression, reuse variable names
    if prefactor
        append!(exprs, [:($(rename(k).name) = cf_from_cl(glq, $(v[1]...), $(rename(k).name); prefactor=true))
                        for (k,v) in cf_table])
    else
        append!(exprs, [:($(rename(k).name) = cf_from_cl(glq, $(v[1]...), $(rename(k).name)))
                        for (k,v) in cf_table])
    end
    # build cf_from_cl expression and add inplace
    for (i,v) in enumerate(cl_table)
        # the awkward pipe is to avoid __dot__ from picking up other variables in the function
        if i == 1
            push!(exprs, :(res = (cl_from_cf(glq,$(v[1]...),lmax, @__dot__$(rename(v[2])))
                                  |> x -> (x .*=@__dot__ $(v[3]); x))))
        else
            push!(exprs, :(res .+= (cl_from_cf(glq,$(v[1]...),lmax, @__dot__$(rename(v[2])))
                                    |> x -> (x .*= @__dot__ $(v[3]); x))))
        end
    end
    exprs
end

function build_l12sum_calculator(expr, name, rename_table, args; prefactor=true, evaluate=false, pre=[], post=[])
    cl_table, cf_table = build_cl_cf_tables(expr; prefactor=prefactor)
    name = name isa String ? Symbol(name) : name
    f = :(function $(name)(lmax, $(map(x->getfield(x,:name), args)...))
              $(pre...)   # allow pass in arbitrary preprocessor
              npoints = (max(lmax,length.([$(args...)])...)*3+1)/2 |> round |> Int
              glq = wignerd.glquad(npoints)
              ℓ = collect(0:(max(length.([$(args...)])...)-1))
              $(build_wigd_calls(cl_table, cf_table, rename_table; prefactor=prefactor)...)
              $(post...)  # allow pass in arbitrary postprocessor
              res         # we have assumed result is stored in this variable
          end)
    # For some reason, functions built this way are contaminated by
    # the types of variables used to built it, i.e., they carry
    # implicite information. Also there is bugs in SymbolicUtils.Code
    # toexpr that it doesn't work with macro like __dot__ nicely.
    # Hence, I use a dirty hack that converts the function to a string
    # and parse it back. This way the implicit type information of
    # all variables are lost in the process.
    evaluate ? (string(linefilter!(f)) |> Meta.parse |> eval) : linefilter!(f)
end
