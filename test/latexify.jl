using Symbolics
using Test
using Latexify

@variables x y z

@test latexify(x^-1) == Latexify.L"$x^{-1}$"

@test latexify((z + x*y^-1) / sin(z)) ==
    Latexify.L"$\left( z + x \cdot y^{-1} \right) \cdot \sin^{-1}\left( z \right)$"

@test latexify((3x - 7y*z^23) * (z - z^2) / x) ==
    Latexify.L"x^{-1} \cdot \left( z -1 \cdot z^{2} \right) \cdot \left( 3 \cdot x -7 \cdot y \cdot z^{23} \right)"
