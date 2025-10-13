module SymbolicsDistributionsExt

using Symbolics
using Distributions
using Symbolics: Num, Arr, VartypeT, unwrap, BasicSymbolic, shape


for f in [pdf, logpdf, cdf, logcdf, quantile]
    @eval function (::$(typeof(f)))(dist::Distributions.Distribution, x::Num)
        $f(dist, unwrap(x))
    end
    @eval function (::$(typeof(f)))(dist::Distributions.Distribution, x::Arr)
        $f(dist, unwrap(x))
    end
    @eval function (::$(typeof(f)))(dist::BasicSymbolic{VartypeT}, x::Num)
        $f(dist, unwrap(x))
    end
    @eval function (::$(typeof(f)))(dist::BasicSymbolic{VartypeT}, x::Arr)
        $f(dist, unwrap(x))
    end
    @eval function (::$(typeof(f)))(dist::BasicSymbolic{VartypeT}, x)
        $f(dist, unwrap(x))
    end
end

for f in [Distributions.Uniform, Distributions.Normal]
    for (T1, T2) in Iterators.product(Iterators.repeated([Real, BasicSymbolic{VartypeT}, Num], 2)...)
        if T1 != Num && T2 != Num
            continue
        end
        @eval function (::Type{$f})(a::$T1, b::$T2)
            $f(unwrap(a), unwrap(b))
        end
    end
end

end
