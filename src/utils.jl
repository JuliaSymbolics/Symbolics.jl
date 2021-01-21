isblock(x) = length(x) == 1 && x[1] isa Expr && x[1].head == :block
function flatten_expr!(x)
    isb = isblock(x)
    if isb
        x = MacroTools.striplines(x[1])
        filter!(z->z isa Symbol || z.head != :line, x.args)
        x = (x.args...,)
    end
    x
end
function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

