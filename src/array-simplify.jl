using Symbolics
import SymbolicUtils: @rule, symtype, Term
import Symbolics: ArrayOp, unwrap, Arr

struct ArrayRule
    r
end

macro array_rule(expr)
    quote
        ArrayRule($(esc(:(@rule($(expr))))))
    end
end

_term(x) = x
_term(x::ArrayOp) = Term{symtype(x.term)}(operation(x.term),
                                          map(_term, arguments(x.term)))

(ar::ArrayRule)(expr::Arr) = wrap(ar(unwrap(expr)))
(ar::ArrayRule)(expr) = ar.r(_term(expr))
