using Symbolics, Aqua
@testset "Aqua" begin
    Aqua.test_ambiguities(Symbolics, recursive = false)
    Aqua.test_deps_compat(Symbolics)
    Aqua.test_piracies(Symbolics,
        treat_as_own = [Symbolics.Symbolic, Symbolics.Sym])
    Aqua.test_project_extras(Symbolics)
    Aqua.test_stale_deps(Symbolics)
    Aqua.test_unbound_args(Symbolics)
    Aqua.test_undefined_exports(Symbolics)
    Aqua.find_persistent_tasks_deps(Symbolics)
end
