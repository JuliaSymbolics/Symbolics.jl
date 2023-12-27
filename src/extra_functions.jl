@register_symbolic Base.binomial(n, k::Integer)::Int true

@register_symbolic Base.sign(x)::Int
derivative(::typeof(sign), args::NTuple{1,Any}, ::Val{1}) = 0

@register_symbolic Base.signbit(x)::Bool
derivative(::typeof(signbit), args::NTuple{1,Any}, ::Val{1}) = 0
derivative(::typeof(abs), args::NTuple{1,Any}, ::Val{1}) = IfElse.ifelse(signbit(args[1]),-one(args[1]),one(args[1]))

function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x < y, one(x), zero(x))
end
function derivative(::typeof(min), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    IfElse.ifelse(x < y, zero(y), one(y))
end
function derivative(::typeof(max), args::NTuple{2,Any}, ::Val{1})
    x, y = args
    IfElse.ifelse(x > y, one(x), zero(x))
end
function derivative(::typeof(max), args::NTuple{2,Any}, ::Val{2})
    x, y = args
    IfElse.ifelse(x > y, zero(y), one(y))
end

@register_symbolic Base.ceil(x)
@register_symbolic Base.floor(x)
@register_symbolic Base.factorial(x)

function derivative(::Union{typeof(ceil),typeof(floor),typeof(factorial)}, args::NTuple{1,Any}, ::Val{1})
    zero(args[1])
end

@register_symbolic Base.rand(x)
@register_symbolic Base.randn(x)

@register_symbolic Distributions.pdf(dist,x)
@register_symbolic Distributions.logpdf(dist,x)
@register_symbolic Distributions.cdf(dist,x)
@register_symbolic Distributions.logcdf(dist,x)
@register_symbolic Distributions.quantile(dist,x)

@register_symbolic Distributions.Uniform(mu,sigma) false
@register_symbolic Distributions.Normal(mu,sigma) false

@register_symbolic ∈(x::Real, y::AbstractArray)::Bool
@register_symbolic ∪(x, y)
@register_symbolic ∩(x, y)
@register_symbolic ∨(x, y)
@register_symbolic ∧(x, y)
@register_symbolic ⊆(x, y)

LinearAlgebra.norm(x::Num, p::Real) = abs(x)
