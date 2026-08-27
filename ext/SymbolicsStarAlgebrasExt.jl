module SymbolicsStarAlgebrasExt

using Symbolics: Differential
import StarAlgebras

Base.:*(D::Differential, x::StarAlgebras.AlgebraElement) = D ∘ x

# StarAlgebras defines scalar multiplication on this non-public abstract type.
Base.:*(D::Differential, x::StarAlgebras.AbstractCoefficients) = D ∘ x

end
