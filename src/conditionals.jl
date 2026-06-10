# `ifelse_eager` and `ifelse_branching` are defined in SymbolicUtils (alongside `ifelse`),
# including their symbolic term construction, type/shape promotion, code generation and
# `simplify` rules. Here we only add the `Num`/`Arr` wrappers, mirroring `ifelse` in `num.jl`
# and `arrays.jl`. Differentiation and linearity/sparsity rules live in `diff.jl`.
for op in (:ifelse_eager, :ifelse_branching)
    @eval begin
        $op(c::Num, x, y) = $op(unwrap(c), unwrap(x), unwrap(y))
        $op(c::Num, x::Num, y) = Num($op(unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x, y::Num) = Num($op(unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Num, y::Num) = Num($op(unwrap(c), unwrap(x), unwrap(y)))

        $op(c::Num, x::Arr{T, N}, y) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Arr{T, N}, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        # Disambiguate mixed scalar/array branches (intersections of the methods above).
        $op(c::Num, x::Num, y::Arr{T, N}) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
        $op(c::Num, x::Arr{T, N}, y::Num) where {T, N} = Arr{T, N}($op(
            unwrap(c), unwrap(x), unwrap(y)))
    end
end
