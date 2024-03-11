struct Struct{T} <: Real
end

"""
    symstruct(T)

Create a symbolic wrapper for struct from a given struct `T`.
"""
Struct(::Type{T}) where T = Struct{T}

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

function symbolic_getproperty(ss, name::Symbol)
    s = symtype(ss)
    idx = findfirst(isequal(name), getelements(s))
    idx === nothing && error("$(juliatype(s)) doesn't have field $(name)!")
    T = getelementtypes(s)[idx]
    SymbolicUtils.term(getfield, ss, Meta.quot(name), type = T)
end
function symbolic_getproperty(s::Union{Arr, Num}, name::Symbol)
    wrap(symbolic_getproperty(unwrap(s), name))
end

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

# We cannot precisely derive the type after `getfield` due to SU limitations,
# so give up and just say Real.
SymbolicUtils.promote_symtype(::typeof(getfield), ::Type{<:Struct}, _...) = Real
SymbolicUtils.promote_symtype(::typeof(setfield!), ::Type{<:Struct}, _, ::Type{T}) where T = T
