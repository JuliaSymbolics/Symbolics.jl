using Symbolics

# https://github.com/JuliaSymbolics/Symbolics.jl/issues/1206
# test there is no extra white spaces on the left of e or f when displaying an array
@testset "e f alignment in array" begin
    r = @variables a b c d e f g h i j

    io = IOBuffer()
    td = TextDisplay(io)

    display(td, r)
    s = String(take!(io))
    @test s == "10-element Vector{Num}:\n a\n b\n c\n d\n e\n f\n g\n h\n i\n j\n"

    m = [a b c d
         b e f g
         c f h i
         d g i j]
    display(td, m)
    s = String(take!(io))
    @test s == "4×4 Matrix{Num}:\n a  b  c  d\n b  e  f  g\n c  f  h  i\n d  g  i  j\n"

    m *= [1, 2, 3, 4]
    display(td, m)
    s = String(take!(io))
    @test s ==
          "4-element Vector{Num}:\n a + 2b + 3c + 4d\n b + 2e + 3f + 4g\n c + 2f + 3h + 4i\n d + 2g + 3i + 4j\n"
end
