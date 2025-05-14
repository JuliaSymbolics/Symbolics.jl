struct Struct{T} <: Real
end

"""
    symstruct(T)

Create a symbolic wrapper for struct from a given struct `T`.
"""
symstruct(::Type{T}) where T = Struct{T}
Struct{T}(vals...) where T = T(vals...)

function Base.hash(x::Struct{T}, seed::UInt) where T
    h1 = hash(T, seed)
    h2 ⊻ (0x0e39036b7de2101a % UInt)
end

"""
    juliatype(s::Type{<:Struct})

Get the Julia type that `s` is representing.
"""
juliatype(::Type{Struct{T}}) where T = T
getelements(s::Type{<:Struct}) = fieldnames(juliatype(s))
getelementtypes(s::Type{<:Struct}) = fieldtypes(juliatype(s))

typed_getfield(obj, ::Val{fieldname}) where fieldname = getfield(obj, fieldname)

"""
    symbolic_getproperty(ss, name::Symbol)

Symbolic term corresponding to accessing the field with name `name`.
"""
function symbolic_getproperty(ss, name::Symbol)
    s = symtype(ss)
    idx = findfirst(isequal(name), getelements(s))
    idx === nothing && error("$(juliatype(s)) doesn't have field $(name)!")
    T = getelementtypes(s)[idx]
    if isstructtype(T)
        T = Struct{T}
    end
    SymbolicUtils.term(typed_getfield, ss, Val{name}(), type = T)
end
function symbolic_getproperty(s::Union{Arr, Num}, name::Symbol)
    wrap(symbolic_getproperty(unwrap(s), name))
end

"""
    symbolic_setproperty!(ss, name::Symbol)

Symbolic term corresponding to modifying the field with name `name` to val `val`.
"""
function symbolic_setproperty!(ss, name::Symbol, val)
    s = symtype(ss)
    idx = findfirst(isequal(name), getelements(s))
    idx === nothing && error("$(juliatype(s)) doesn't have field $(name)!")
    T = getelementtypes(s)[idx]
    SymbolicUtils.term(setfield!, ss, Meta.quot(name), val, type = T)
end
function symbolic_setproperty!(s::Union{Arr, Num}, name::Symbol, val)
    wrap(symbolic_setproperty!(unwrap(s), name, val))
end

function symbolic_constructor(s::Type{<:Struct}, vals...)
    N = length(getelements(s))
    N′ = length(vals)
    N′ == N || error("$(juliatype(s)) needs $N field. Got $N′ fields!")
    SymbolicUtils.term(s, vals..., type = s)
end

# We cannot precisely derive the type after `getfield` due to SU limitations,
# so give up and just say Real.
function SymbolicUtils.promote_symtype(::typeof(typed_getfield), ::Type{<:Struct{T}}, v::Type{Val{fieldname}}) where {T, fieldname} 
    FT = fieldtype(T, fieldname)
    if isstructtype(FT)
        return Struct{FT}
    end 
    FT
end

function SymbolicUtils.promote_symtype(s::Type{<:Struct{T}}, _...) where T
    s
end

SymbolicUtils.promote_symtype(::typeof(setfield!), ::Type{<:Struct}, _, ::Type{T}) where T = T
