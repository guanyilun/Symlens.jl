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

reverse(x::Dict) =  Dict(v=>k for (k,v) in x)

function is_atomic(expr)
    vars = get_variables_deep(expr)
    length(vars)>1 ? false : true
end

function are_factors_atomic(expr)
    all(map(ex->is_atomic(ex), arguments(expr)))
end

is_term_of(expr, syms) = (vars = get_variables_deep(expr); Set(vars) == Set(syms))

count_terms_of(exprs::Vector, target::SymbolicUtils.Sym{SymbolicUtils.FnType{T,N},Nothing}) where {T,N} = (x-> x isa SymbolicUtils.Term ? x.f==target : false).(exprs) |> sum
