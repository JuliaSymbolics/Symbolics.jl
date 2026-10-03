# Linear expansion

Optimizations of matrix `linear_expansion` must preserve the variable-ordered `LinearExpander` behavior. Later unknowns can appear in earlier coefficient columns, and repeated or overlapping unknowns are order-sensitive. Use the original expansion path whenever a fast path cannot establish equivalent results, including products with cancelling factors.

Cover both expansion and `symbolic_linear_solve` in `test/linear_solver.jl`, and compare against the unmodified base for composite unknowns, cancellation, and array expressions. Run the full file and `GROUP=Core julia --project=. -e 'using Pkg; Pkg.test()'` before pushing solver changes.
