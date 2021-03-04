using ModelingToolkit, ParameterizedFunctions
f = @ode_def_bare LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a b c d
p = [1.5,1,3,1]
tspan = (0.0,10.0)
u0 = [1.0,1.0]
prob = ODEProblem(f,u0,tspan,p)

sys = modelingtoolkitize(prob)
[eq.rhs for eq ∈ equations(sys)]
matstr = Symbolics.build_function(map(x->x.rhs,equations(sys)),states(sys),
                                        parameters(sys),independent_variable(sys),
                                        target = ModelingToolkit.MATLABTarget())
@test matstr == "diffeqf = @(t,internal_var___u) [
  internal_var___p(1) * internal_var___u(1) - internal_var___p(2) * internal_var___u(1) * internal_var___u(2);
  internal_var___u(2) * -(internal_var___p(3)) + internal_var___p(4) * internal_var___u(1) * internal_var___u(2);
];
"
