using Symbolics
import Symbolics: value, coeff, sympy_integrate

function laplace(f, t, s)
    # from https://tutorial.math.lamar.edu/Classes/DE/Laplace_Table.aspx
    transform_rules = [
        @rule 1 => 1/s
        @rule exp(t) => 1/(s - 1)
        @rule exp(~a * t) => 1/(s - ~a)
        @rule t => 1/s^2
        @rule t^~n => factorial(~n)/s^(~n + 1)
        @rule sqrt(t) => sqrt(pi)/(2 * s^(3/2))
        @rule sin(t) => 1/(s^2 + 1)
        @rule sin(~a * t) => ~a/(s^2 + ~a^2)
        @rule cos(t) => s/(s^2 + 1)
        @rule cos(~a * t) => s/(s^2 + ~a^2)
        @rule t*sin(t) => 1/(s^2 + 1)^2
        @rule t*sin(~a * t) => 2*~a*s / (s^2 + ~a^2)^2
        @rule t*cos(t) => (s^2 - 1) / (s^2 + 1)^2
        @rule t*cos(~a * t) => (s^2 - ~a^2) / (s^2 + ~a^2)^2
        @rule sin(t) - t*cos(t) => 2 / (s^2 + 1)^2
        @rule sin(~a*t) - ~a*t*cos(~a*t) => 2*~a^3 / (s^2 + ~a^2)^2
        @rule sin(t) + t*cos(t) => 2s^2 / (s^2 + 1)^2
        @rule sin(~a*t) + ~a*t*cos(~a*t) => 2*~a*s^2 / (s^2 + ~a^2)^2
        @rule cos(~a*t) - ~a*t*sin(~a*t) => s*(s^2 + ~a^2) / (s^2 + ~a^2)^2
        @rule cos(~a*t) + ~a*t*sin(~a*t) => s*(s^2 + 3*~a^2) / (s^2 + ~a^2)^2
        @rule sin(~a*t + ~b) => (s*sin(~b) + ~a*cos(~b)) / (s^2 + ~a^2)
        @rule cos(~a*t + ~b) => (s*cos(~b) - ~a*sin(~b)) / (s^2 + ~a^2)
        @rule sinh(~a * t) => ~a/(s^2 - ~a^2)
        @rule cosh(~a * t) => s/(s^2 - ~a^2)
        @rule exp(~a*t) * sin(~b * t) => ~b / ((s-~a)^2 + ~b^2)
        @rule exp(~a*t) * cos(~b * t) => (s-~a) / ((s-~a)^2 + ~b^2)
        @rule exp(~a*t) * sinh(~b * t) => ~b / ((s-~a)^2 - ~b^2)
        @rule exp(~a*t) * cosh(~b * t) => (s-~a) / ((s-~a)^2 - ~b^2)
        @rule t^~n * exp(~a * t) => factorial(~n) / (s - ~a)^(~n + 1)
        @rule exp(~c*t) * ~g => laplace(~g, t, s - ~c)
    ]



    return sympy_integrate(f * exp(-s * t), (t, 0, Inf))
end

function inverse_laplace(F, t, s)
    inverse_transform_rules = [
        @rule 1/s => 1
        @rule 1/(s + ~a) => exp(~a * t)
        @rule factorial(~n)/s^(~n + 1) => t^~n
        @rule sqrt(pi)/(2 * s^(3/2)) => sqrt(t)
        @rule ~a/(s^2 + ~a^2) => sin(~a * t)
        @rule s/(s^2 + ~a^2) => cos(~a * t)
        @rule 2*~a*s / (s^2 + ~a^2)^2 => t*sin(~a * t)
        @rule (s^2 - ~a^2) / (s^2 + ~a^2)^2 => t*cos(~a * t)
        @rule 2*~a^3 / (s^2 + ~a^2)^2 => sin(~a*t) - ~a*t*cos(~a*t)
        @rule 2*~a*s^2 / (s^2 + ~a^2)^2 => sin(~a*t) + ~a*t*cos(~a*t)
        @rule s*(s^2 + ~a^2) / (s^2 + ~a^2)^2 => cos(~a*t) - ~a*t*sin(~a*t)
        @rule s*(s^2 + 3*~a^2) / (s^2 + ~a^2)^2 => cos(~a*t) + ~a*t*sin(~a*t)
        @rule (s*sin(~b) + ~a*cos(~b)) / (s^2 + ~a^2) => sin(~a*t + ~b)
        @rule (s*cos(~b) - ~a*sin(~b)) / (s^2 + ~a^2) => cos(~a*t + ~b)
        @rule ~b/(s^2 - ~b^2) => sinh(~b * t)
        @rule s/(s^2 - ~b^2) => cosh(~b * t)
        @rule ~b / ((s-~c)^2 + ~b^2) => exp(~c*t) * sin(~b * t)
        @rule (s-~c) / ((s-~c)^2 + ~b^2) => exp(~c*t) * cos(~b * t)
        @rule ~b / ((s-~c)^2 - ~b^2) => exp(~c*t) * sinh(~b * t)
        @rule (s-~c) / ((s-~c)^2 - ~b^2) => exp(~c*t) * cosh(~b * t)
        @rule factorial(~n) / (s - ~a)^(~n + 1) => t^~n * exp(~a * t)
    ]
end