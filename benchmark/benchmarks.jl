using BenchmarkTools, Symbolics

using Random, LinearAlgebra

SUITE = BenchmarkGroup()

vars = @variables a,b,c,d,e,f,g,h,i


X = [0 b c;
     d e f;
     g h i]

z = det(X) - (-b * (d*i-f*g) + c * (d*h - e*g))

F = lu(X)
FbyX = F \ X

x = (f + ((((g*(c^2)*(e^2)) / d - e*h*(c^2)) / b + (-c*e*f*g) / d + c*e*i) /
          (i + ((c*e*g) / d - c*h) / b + (-f*g) / d) - c*e) / b +
     ((g*(f^2)) / d + ((-c*e*f*g) / d + c*f*h) / b - f*i) /
     (i + ((c*e*g) / d - c*h) / b + (-f*g) / d)) / d

o = (d + (e*((c*(g + (-d*g) / d)) / (i + (-c*(h + (-e*g) / d)) / b + (-f*g) / d))) / b + (-f*(g + (-d*g) / d)) / (i + (-c*(h + (-e*g) / d)) / b + (-f*g) / d)) / d

SUITE["iszero/1"] = @benchmarkable iszero($((b*(h + (-e*g) / d)) / b + (e*g) / d - h))
SUITE["isone/1"] = @benchmarkable FbyX == I

SUITE["iszero/2"] = @benchmarkable iszero(x)
SUITE["isone/2"] = @benchmarkable isone(o)
SUITE["lu"] = @benchmarkable lu(X)

SUITE["_solve"] = @benchmarkable Symbolics._solve(X, X, true)
