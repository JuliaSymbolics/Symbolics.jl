# Transfer Functions and Frequency Response

The transfer function of a linear state-space system `x' = Ax + Bu`, `y = Cx + Du` can be
computed symbolically, keeping the parameters, and then evaluated on a frequency grid. For
example, for a series RLC circuit with the capacitor voltage and the inductor current as
states, driven by a voltage source and measured across the capacitor:

```@example transfer
using Symbolics
@variables s R L C
A = [0 1/C; -1/L -R/L]
B = [0, 1/L]
G = transfer_function(A, B, [1, 0], 0, s)
```

Substitute values for the parameters to obtain the frequency response at `s = im*ω`:

```@example transfer
H = frequency_response(substitute(G, Dict(R => 1.0, L => 0.5, C => 0.2)), s, [0.1, 1.0, 10.0])
```

Magnitude and phase of the first input-output pair are `abs.(H[1, 1, :])` and
`angle.(H[1, 1, :])`.

The matrices of a nonlinear model linearized symbolically, for example by
ModelingToolkit's `linearize_symbolic`, can be passed directly. That linearization uses
symbolic differentiation, so it only works for models whose operations can all be
differentiated symbolically. Systems in descriptor form `Ex' = Ax + Bu` are not supported;
pass the reduced matrices.

```@docs
transfer_function
frequency_response
```
