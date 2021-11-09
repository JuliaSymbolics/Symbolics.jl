# optimize(o) |> Symbolics.unwrap |> x -> Postwalk(y -> istree(y) ? similarterm(SymbolicUtils.Symbolic, operation(y), arguments(y), symtype(y); metadata=metadata(y)) : y)(x)

using BenchmarkTools, SymbolicUtils
using Symbolics
using Symbolics: optimize, build_function, get_variables

import SymbolicUtils.Code: toexpr

ops = [*, +, -]
syms = (@variables a b c d e f g h i j k l m ) |> collect
# n o p q r s t u v w x y z

constant_vals = -10:0.3:10
function rand_expr(d=0, mindepth=70, maxdepth=100, ops=ops)
    if d < mindepth
        (rand(ops))(rand_expr(d+1, mindepth, maxdepth), rand_expr(d+1, mindepth, maxdepth))
    elseif d >= maxdepth
        rand(0:8) == 0 ?  rand(constant_vals) : rand([-1,1]) * rand(syms)
    elseif mindepth <= d < maxdepth
        rand(Bool) ? (rand(0:2) == 0 ? rand(constant_vals) : rand(syms)) : term(rand(ops), rand_expr(d+1, mindepth, maxdepth), rand_expr(d+1, mindepth, maxdepth))
    end 
end


# x = (f + ((((g*(c^2)*(e^2)) / d - e*h*(c^2)) / b + (-c*e*f*g) / d + c*e*i) /
#           (i + ((c*e*g) / d - c*h) / b + (-f*g) / d) - c*e) / b +
#      ((g*(f^2)) / d + ((-c*e*f*g) / d + c*f*h) / b - f*i) /
#      (i + ((c*e*g) / d - c*h) / b + (-f*g) / d)) / d


# getcost(x)

# x_expr = toexpr(x)
# opt = optimize(x)
# getcost(opt)
# opt_expr = toexpr(opt)
tofun(ex;kws...) = build_function(ex, syms...; kws...) |> eval

x = rand_expr(0, 12, 12)



# @btime ($(tofun(x)))($v...)
# @btime ($(tofun(x; cycle_optimize=true)))($v...)
f = tofun(x; cycle_optimize=false) 
f_opt = tofun(x; cycle_optimize=true) 

v = rand(13)
@btime ($f)($v...)
@btime ($f_opt)($v...)

# @code_lowered f(v...)

using Suppressor


a = 0
b = 0
c = 0
d = 0

for i = 1:10
    x = rand_expr(0, 12, 12)

    # @btime ($(tofun(x)))($v...)
    # @btime ($(tofun(x; cycle_optimize=true)))($v...)
    f = tofun(x; cycle_optimize=false) 
    f_opt = tofun(x; cycle_optimize=true) 

    v = rand(13)

    llvm_code_unopt = @capture_out begin 
        @code_llvm optimize=false debuginfo=:none f(v...)
    end;
    a += length(split(llvm_code_unopt, "\n"))

    llvm_code_unopt_eqsat = @capture_out begin 
        @code_llvm optimize=false debuginfo=:none f_opt(v...)
    end;
    b += length(split(llvm_code_unopt_eqsat, "\n"))

    llvm_code_opt = @capture_out begin 
        @code_llvm debuginfo=:none f(v...)
    end;
    c += length(split(llvm_code_opt, "\n"))

    llvm_code_opt_eqsat = @capture_out begin 
        @code_llvm debuginfo=:none f_opt(v...)
    end;
    d += length(split(llvm_code_opt_eqsat, "\n"))
end
# clipboard(@code_lowered f_opt(v...))
a/10 
b/10 
c/10
d/10


ops = [*, +, -]
