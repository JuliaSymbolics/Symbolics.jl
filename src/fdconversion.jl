module FDConversion

import Symbolics
import SymbolicUtils
import FastDifferentiation as FD
import Random

#=
try implementing these
    https://github.com/JuliaSymbolics/SymbolicUtils.jl/blob/master/src/interface.jl
=#

"""converts from Node to Symbolics expression"""
function _to_symbolics!(a::T, cache::IdDict, variable_map::IdDict) where {T<:FD.Node}
    tmp = get(cache, a, nothing)
    if tmp !== nothing
        return tmp
    else
        if FD.arity(a) === 0
            if FD.is_constant(a)
                cache[a] = Symbolics.Num(FD.value(a))
            elseif FD.is_variable(a)
                tmp = Symbolics.variable(FD.value(a))
                cache[a] = tmp
                variable_map[a] = tmp
            else
                throw(ErrorException("Node with 0 children was neither constant nor variable. This should never happen."))
            end
        else
            if FD.arity(a) === 1
                cache[a] = a.node_value(_to_symbolics!(a.children[1], cache, variable_map))
            else
                cache[a] = foldl(a.node_value, _to_symbolics!.(a.children, Ref(cache), Ref(variable_map)))
            end
        end
    end
end

function to_symbolics(a::T, cache::IdDict=IdDict(), variable_map::IdDict=IdDict()) where {T<:FD.Node}
    _to_symbolics!(a, cache, variable_map)
    return cache[a], variable_map
end
export to_symbolics

function to_symbolics(a::AbstractArray{T}, cache::IdDict=IdDict(), variable_map::IdDict=IdDict()) where {T<:FD.Node}
    _to_symbolics!.(a, Ref(cache), Ref(variable_map))
    return map(x -> cache[x], a), variable_map
end

"""Converts expression x from Symbolics to FastDifferentiation form. Returns a tuple of the **FD** form of the expression and a Dict that maps Symbolics variables to **FD** variables"""
function to_fd(x::Real, cache=IdDict(), substitutions=IdDict())
    result = _to_FD!(x, cache, substitutions)
    syms = collect(filter(x -> SymbolicUtils.issym(x), keys(cache)))
    sym_map = Dict(zip(syms, map(x -> cache[x], syms)))

    return result, sym_map
end
export to_fd

function to_fd(x::AbstractArray{<:Real})
    cache = IdDict()
    substitutions = IdDict()
    result = _to_FD!.(x, Ref(cache), Ref(substitutions))
    syms = collect(filter(x -> SymbolicUtils.issym(x), keys(cache)))
    sym_map = Dict(zip(syms, map(x -> cache[x], syms))) #map from symbolics variables to FD variables
    return result, sym_map
end

function _to_FD!(sym_node, cache::IdDict, visited::IdDict)
    # Substitutions are done on a Node graph, not a SymbolicsUtils.Sym graph. As a consequence the values
    # in substitutions Dict are Node not Sym type. cache has keys (op,args...) where op is generally a function type but sometimes a Sym, 
    # and args are all Node types.

    @assert !(typeof(sym_node) <: Symbolics.Arr) "Differentiation of expressions involving arrays and array variables is not yet supported."
    if SymbolicUtils.istree(Symbolics.unwrap(sym_node))
        @assert !SymbolicUtils.issym(SymbolicUtils.operation(Symbolics.unwrap(sym_node))) "expressions of the form `q(t)` not yet supported."
    end

    symx = isa(sym_node, Symbolics.Num) ? sym_node.val : sym_node
    @assert typeof(symx) != Symbolics.Num


    tmpsub = get(visited, symx, nothing)
    if tmpsub !== nothing
        return visited[symx] #substitute Node object for symbolic object
    end

    tmp = get(cache, symx, nothing)

    if tmp !== nothing
        return tmp
    elseif !SymbolicUtils.istree(symx)
        if SymbolicUtils.issym(symx)
            tmpnode = FD.Node(Symbol(symx))
        else #must be a number of some kind
            tmpnode = FD.Node(symx)
        end

        cache[symx] = tmpnode

        return tmpnode
    else
        numargs = length(SymbolicUtils.arguments(symx))
        symargs = SymbolicUtils.arguments(symx)

        args = _to_FD!.(symargs, Ref(cache), Ref(visited))

        key = (SymbolicUtils.operation(symx), args...)
        tmp = get(cache, key, nothing)

        if tmp !== nothing
            return tmp
        else
            tmpnode = FD.Node(SymbolicUtils.operation(symx), args...)
            cache[key] = tmpnode

            return tmpnode
        end
    end
end


function remap(fd_function_to_call, symbolics_function, differentiation_variables::AbstractVector{Symbolics.Num}, fast_differentiation::Bool)
    fd_func, variable_map = to_fd(symbolics_function)
    partial_vars = map(x -> variable_map[x], differentiation_variables)
    tmp = fd_function_to_call(fd_func, partial_vars)
    if fast_differentiation
        return fd_func, variable_map #return FastDifferentiation expression to be passed to make_function for efficient evaluation
    else
        reverse_map = IdDict{Any,Any}(map(x -> variable_map[x] => x, differentiation_variables))
        reverse_variable_map = deepcopy(reverse_map)
        return to_symbolics(tmp, reverse_map, reverse_variable_map) #return Symbolics expression for further evaluation
    end
end


"""
Converts from `Symbolics` form to `FastDifferentiation` form and computes Jacobian with respect to `diff_variables`.
If `fast_differentiation=false` returns result in Symbolics form. If `fast_differentiation=true` the result will be a 2-tuple. The first tuple entry will be the jacobian of `symbolics_function` converted to `FastDifferentiation` form. 
The second tuple term will be `differentiation_variables` converted to `FastDifferentiation` form.
This tuple can be passed to `fd_make_function` to create an efficient executable."""
fd_jacobian(symbolics_function::AbstractArray{Symbolics.Num}, differentiation_variables::AbstractVector{Symbolics.Num}; fast_differentiation=false) = remap(FD.jacobian, symbolics_function, differentiation_variables, fast_differentiation)
export fd_jacobian

"""
Converts from `Symbolics` form to `FastDifferentiation` form and computes sparse Jacobian with respect to `diff_variables`.
If `fast_differentiation=false` returns result in Symbolics form. If `fast_differentiation=true` the result will be a 2-tuple. The first tuple entry will be the sparse jacobian of `symbolics_function` converted to `FastDifferentiation` form. 
The second tuple term will be `differentiation_variables` converted to `FastDifferentiation` form.
This tuple can be passed to `fd_make_function` to create an efficient executable."""
fd_sparse_jacobian(symbolics_function::AbstractArray{Symbolics.Num}, differentiation_variables::AbstractVector{Symbolics.Num}; fast_differentiation=false) = remap(FD.sparse_jacobian, symbolics_function, differentiation_variables, fast_differentiation)
export fd_sparse_jacobian

"""
Converts from `Symbolics` form to `FastDifferentiation` form and computes Hessian with respect to `diff_variables`.
If `fast_differentiation=false` returns result in Symbolics form. If `fast_differentiation=true` the result will be a 2-tuple. The first tuple entry will be the hessian of `symbolics_function` converted to `FastDifferentiation` form. 
The second tuple term will be `differentiation_variables` converted to `FastDifferentiation` form.
This tuple can be passed to `fd_make_function` to create an efficient executable."""
fd_hessian(symbolics_function::Symbolics.Num, differentiation_variables::AbstractVector{Symbolics.Num}; fast_differentiation=false) = remap(FD.hessian, symbolics_function, differentiation_variables, fast_differentiation)
export fd_sparse_hessian

"""
Converts from `Symbolics` form to `FastDifferentiation` form and computes sparse Hessian with respect to `diff_variables`.
If `fast_differentiation=false` returns result in Symbolics form. If `fast_differentiation=true` the result will be a 2-tuple. The first tuple entry will be the sparse Hessian of `symbolics_function` converted to `FastDifferentiation` form. 
The second tuple term will be `differentiation_variables` converted to `FastDifferentiation` form.
This tuple can be passed to `fd_make_function` to create an efficient executable."""
fd_sparse_hessian(symbolics_function::Symbolics.Num, differentiation_variables::AbstractVector{Symbolics.Num}; fast_differentiation=false) = remap(FD.sparse_hessian, symbolics_function, differentiation_variables, fast_differentiation)
export fd_sparse_hessian

"""
Creates an efficient runtime generated function to evaluate the function defined in `a[1]` with arguments given by `a[2]`. Used in conjuction with fd_jacobian,fd_sparse_jacobian,fd_hessian,fd_sparse_hessian.
"""
fd_make_function(a::Tuple{T,S}; in_place=false) where {T<:AbstractArray{<:FD.Node},S<:AbstractVector{<:FD.Node}} = FD.make_function(a[1], a[2], in_place=in_place)
export fd_make_function
#etc. for Jv Jᵀv Hv

# export FastDifferentiation.make_function

end # module FSDConvert

