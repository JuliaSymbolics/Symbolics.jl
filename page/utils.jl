# Thanks @tlienart

using Markdown

function hfun_doc(params)
    fname = params[1]
    head = length(params) > 1 ? params[2] : fname
    type = length(params) == 3 ? params[3] : ""
    mod = length(params) == 4 ? params[4] : ""
    if mod != ""
        doc = eval(Meta.parse("using SymbolicUtils; @doc SymbolicUtils.$mod.$fname"))
    else
        doc = eval(Meta.parse("using SymbolicUtils; @doc SymbolicUtils.$fname"))
    end
    txt = Markdown.plain(doc)
    # jldoctest blocks don't get syntax highlighting in Franklin.jl.
    txt = replace(txt, "```jldoctest" => "```")
    # possibly further processing here
    body = Franklin.fd2html(txt, internal=true)
    return """
      <div class="docstring">
          <h2 class="doc-header" id="$fname">
            <a href="#$fname">$head</a>
            <div class="doc-type">$type</div></h2>
          <div class="doc-content">$body</div>
      </div>
    """
end

function stringmime(m, x)
    open("/tmp/foo", "w") do io
        show(io, m, x)
    end
    sprint(io->show(io, m, x))
end

const mime_renderers = [
    MIME("text/latex") => (x,m) -> stringmime(m, x),
    MIME("text/markdown") => (x,m) -> stringmime(m, x),
    MIME("text/html") => (x,m) -> "~~~$(stringmime(m, x))~~~",
    MIME("image/svg+xml") => (x,m) -> "~~~$(stringmime(m, x))~~~",
    MIME("image/png") => (x,m) -> "cant render png yet",
    MIME("image/jpeg") => (x,m) -> "cant render jpeg yet",
    MIME("text/plain") => (x,m) -> "`$x`",
]

function print_bestmime(x)
    if x === nothing
        return print()
    end
    for (m, f) in mime_renderers
        if showable(m, x)
            print(f(x, m))
            return
        end
    end
    print("`could not render the result`")
end

function repl_cell(ex, hide_output)
    s = """```julia:repl-cell
    res = begin # hide
    $ex
    end # hide
    Franklin.utils_module().print_bestmime(res) # hide
    ```
    """
    if !hide_output
        s *= """~~~<div class="cell-output">~~~\\textoutput{repl-cell}~~~</div>~~~"""
    end
    s
end

function lx_repl(com, _) # the signature must look like this
    # leave this first line, it extracts the content of the brace
    content = Franklin.content(com.braces[1])
    # dumb way to recover stuff
    lines = split(content, "\n")
    cells = join([repl_cell(l, endswith(strip(l), ";"))  for l in lines if !isempty(strip(l)) ], "\n")
    str = """~~~<div class="repl-block">~~~$cells~~~</div>~~~"""
    print(str)
    str
end
