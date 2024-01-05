module SymbolicsForwardDiffExt

using ForwardDiff
using ForwardDiff.NaNMath
using ForwardDiff.DiffRules
using ForwardDiff: value, Dual, partials
using Symbolics

const AMBIGUOUS_TYPES = (Num,)

####################################
# N-ary Operation Definition Tools #
####################################

macro define_binary_dual_op(f, xy_body, x_body, y_body, Ts)
    FD = ForwardDiff
    defs = quote end
    for R in Ts
        expr = quote
            @inline $(f)(x::$FD.Dual{Tx}, y::$R) where {Tx} = $x_body
            @inline $(f)(x::$R, y::$FD.Dual{Ty}) where {Ty} = $y_body
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

macro define_ternary_dual_op(f, xyz_body, xy_body, xz_body, yz_body, x_body, y_body, z_body, Ts)
    FD = ForwardDiff
    defs = quote end
    for R in Ts
        expr = quote
            @inline $(f)(x::$FD.Dual{Txy}, y::$FD.Dual{Txy}, z::$R) where {Txy} = $xy_body
            @inline $(f)(x::$FD.Dual{Tx}, y::$FD.Dual{Ty}, z::$R)  where {Tx, Ty} = Ty ≺ Tx ? $x_body : $y_body
            @inline $(f)(x::$FD.Dual{Txz}, y::$R, z::$FD.Dual{Txz}) where {Txz} = $xz_body
            @inline $(f)(x::$FD.Dual{Tx}, y::$R, z::$FD.Dual{Tz}) where {Tx,Tz} = Tz ≺ Tx ? $x_body : $z_body
            @inline $(f)(x::$R, y::$FD.Dual{Tyz}, z::$FD.Dual{Tyz}) where {Tyz} = $yz_body
            @inline $(f)(x::$R, y::$FD.Dual{Ty}, z::$FD.Dual{Tz}) where {Ty,Tz} = Tz ≺ Ty ? $y_body : $z_body
        end
        append!(defs.args, expr.args)
        for Q in Ts
            Q === R && continue
            expr = quote
                @inline $(f)(x::$FD.Dual{Tx}, y::$R, z::$Q) where {Tx} = $x_body
                @inline $(f)(x::$R, y::$FD.Dual{Ty}, z::$Q) where {Ty} = $y_body
                @inline $(f)(x::$R, y::$Q, z::$FD.Dual{Tz}) where {Tz} = $z_body
            end
            append!(defs.args, expr.args)
        end
        expr = quote
            @inline $(f)(x::$FD.Dual{Tx}, y::$R, z::$R) where {Tx} = $x_body
            @inline $(f)(x::$R, y::$FD.Dual{Ty}, z::$R) where {Ty} = $y_body
            @inline $(f)(x::$R, y::$R, z::$FD.Dual{Tz}) where {Tz} = $z_body
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

function binary_dual_definition(M, f, Ts)
    FD = ForwardDiff
    dvx, dvy = DiffRules.diffrule(M, f, :vx, :vy)
    Mf = M == :Base ? f : :($M.$f)
    xy_work = FD.qualified_cse!(quote
        val = $Mf(vx, vy)
        dvx = $dvx
        dvy = $dvy
    end)
    dvx, _ = DiffRules.diffrule(M, f, :vx, :y)
    x_work = FD.qualified_cse!(quote
        val = $Mf(vx, y)
        dvx = $dvx
    end)
    _, dvy = DiffRules.diffrule(M, f, :x, :vy)
    y_work = FD.qualified_cse!(quote
        val = $Mf(x, vy)
        dvy = $dvy
    end)
    expr = quote
        @define_binary_dual_op(
            $M.$f,
            begin
                vx, vy = $FD.value(x), $FD.value(y)
                $xy_work
                return $FD.dual_definition_retval(Val{Txy}(), val, dvx, $FD.partials(x), dvy, $FD.partials(y))
            end,
            begin
                vx = $FD.value(x)
                $x_work
                return $FD.dual_definition_retval(Val{Tx}(), val, dvx, $FD.partials(x))
            end,
            begin
                vy = $FD.value(y)
                $y_work
                return $FD.dual_definition_retval(Val{Ty}(), val, dvy, $FD.partials(y))
            end,
            $Ts
        )
    end
    return expr
end

###################################
# General Mathematical Operations #
###################################

for (M, f, arity) in DiffRules.diffrules(filter_modules = nothing)
    if (M, f) in ((:Base, :^), (:NaNMath, :pow), (:Base, :/), (:Base, :+), (:Base, :-), (:Base, :sin), (:Base, :cos))
        continue  # Skip methods which we define elsewhere.
    elseif !(isdefined(@__MODULE__, M) && isdefined(getfield(@__MODULE__, M), f))
        continue  # Skip rules for methods not defined in the current scope
    end
    if arity == 1
        # no-op
    elseif arity == 2
        eval(binary_dual_definition(M, f, AMBIGUOUS_TYPES))
    else
        # error("ForwardDiff currently only knows how to autogenerate Dual definitions for unary and binary functions.")
        # However, the presence of N-ary rules need not cause any problems here, they can simply be ignored.
    end
end

#################
# Special Cases #
#################

# +/- #
#-----#

@eval begin
    @define_binary_dual_op(
        Base.:+,
        begin
            vx, vy = value(x), value(y)
            Dual{Txy}(vx + vy, partials(x) + partials(y))
        end,
        Dual{Tx}(value(x) + y, partials(x)),
        Dual{Ty}(x + value(y), partials(y)),
        $AMBIGUOUS_TYPES
    )
end

@eval begin
    @define_binary_dual_op(
        Base.:-,
        begin
            vx, vy = value(x), value(y)
            Dual{Txy}(vx - vy, partials(x) - partials(y))
        end,
        Dual{Tx}(value(x) - y, partials(x)),
        Dual{Ty}(x - value(y), -partials(y)),
        $AMBIGUOUS_TYPES
    )
end

# / #
#---#

# We can't use the normal diffrule autogeneration for this because (x/y) === (x * (1/y))
# doesn't generally hold true for floating point; see issue #264
@eval begin
    @define_binary_dual_op(
        Base.:/,
        begin
            vx, vy = value(x), value(y)
            Dual{Txy}(vx / vy, _div_partials(partials(x), partials(y), vx, vy))
        end,
        Dual{Tx}(value(x) / y, partials(x) / y),
        begin
            v = value(y)
            divv = x / v
            Dual{Ty}(divv, -(divv / v) * partials(y))
        end,
        $AMBIGUOUS_TYPES
    )
end

# exponentiation #
#----------------#

for f in (:(Base.:^), :(NaNMath.pow))
    @eval begin
        @define_binary_dual_op(
            $f,
            begin
                vx, vy = value(x), value(y)
                expv = ($f)(vx, vy)
                powval = vy * ($f)(vx, vy - 1)
                if isconstant(y)
                    logval = one(expv)
                elseif iszero(vx) && vy > 0
                    logval = zero(vx)
                else
                    logval = expv * log(vx)
                end
                new_partials = _mul_partials(partials(x), partials(y), powval, logval)
                return Dual{Txy}(expv, new_partials)
            end,
            begin
                v = value(x)
                expv = ($f)(v, y)
                if y == zero(y) || iszero(partials(x))
                    new_partials = zero(partials(x))
                else
                    new_partials = partials(x) * y * ($f)(v, y - 1)
                end
                return Dual{Tx}(expv, new_partials)
            end,
            begin
                v = value(y)
                expv = ($f)(x, v)
                deriv = (iszero(x) && v > 0) ? zero(expv) : expv*log(x)
                return Dual{Ty}(expv, deriv * partials(y))
            end,
            $AMBIGUOUS_TYPES
        )
    end
end

# hypot #
#-------#

@eval begin
    @define_ternary_dual_op(
        Base.hypot,
        calc_hypot(x, y, z, Txyz),
        calc_hypot(x, y, z, Txy),
        calc_hypot(x, y, z, Txz),
        calc_hypot(x, y, z, Tyz),
        calc_hypot(x, y, z, Tx),
        calc_hypot(x, y, z, Ty),
        calc_hypot(x, y, z, Tz),
        $AMBIGUOUS_TYPES
    )
end

# fma #
#-----#

@eval begin
    @define_ternary_dual_op(
        Base.fma,
        calc_fma_xyz(x, y, z),                         # xyz_body
        calc_fma_xy(x, y, z),                          # xy_body
        calc_fma_xz(x, y, z),                          # xz_body
        Base.fma(y, x, z),                             # yz_body
        Dual{Tx}(fma(value(x), y, z), partials(x) * y), # x_body
        Base.fma(y, x, z),                              # y_body
        Dual{Tz}(fma(x, y, value(z)), partials(z)),     # z_body
        $AMBIGUOUS_TYPES
    )
end

# muladd #
#--------#

@eval begin
    @define_ternary_dual_op(
        Base.muladd,
        calc_muladd_xyz(x, y, z),                         # xyz_body
        calc_muladd_xy(x, y, z),                          # xy_body
        calc_muladd_xz(x, y, z),                          # xz_body
        Base.muladd(y, x, z),                             # yz_body
        Dual{Tx}(muladd(value(x), y, z), partials(x) * y), # x_body
        Base.muladd(y, x, z),                             # y_body
        Dual{Tz}(muladd(x, y, value(z)), partials(z)),     # z_body
        $AMBIGUOUS_TYPES
    )
end

end
