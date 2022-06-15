using Symbolics, Test, SparseArrays
@variables a b c

# Auxiliary Functions
get_sparsity_pattern(h::Union{SparseVector{Num}, SparseMatrixCSC{Num,Int}}) = get_sparsity_pattern(Array(h))
get_sparsity_pattern(h::Array{Num}) = sparse(Int.(.!isequal.(h, 0)))

input = [1, 2, 3]

# ===== Dense tests =====
# Arrays of Matrices
h_dense_arraymat = [[a 1; b 0], [0 0; 0 0], [a c; 1 0]] # empty array support required
function h_dense_arraymat_julia!(out, x)
    a, b, c = x
    out[1] .= [a[1] 1; b[1] 0]
    out[2] .= [0 0; 0 0]
    out[3] .= [a[1] c[1]; 1 0]
end

h_dense_arraymat_ip! = eval(Symbolics.build_function(h_dense_arraymat, [a, b, c])[2])
h_dense_arraymat_ip_skip! = eval(Symbolics.build_function(h_dense_arraymat, [a, b, c], skipzeros=true, fillzeros=false)[2])
out_1_arraymat = [fill(42, 2, 2) for i in 1:3]
out_2_arraymat = deepcopy(out_1_arraymat)
h_dense_arraymat_julia!(out_1_arraymat, input)
h_dense_arraymat_ip_skip!(out_2_arraymat, input)
@test all(isequal(42), out_2_arraymat[2])
foreach(mat->fill!(mat, 0), out_2_arraymat)
h_dense_arraymat_ip_skip!(out_2_arraymat, input)
@test out_1_arraymat == out_2_arraymat
h_dense_arraymat_ip!(out_2_arraymat, input)
@test out_1_arraymat == out_2_arraymat

# Arrays of Matrices, Heterogeneous element types
test_exp = [exp(a) * exp(b), a]
h_dense_arraymat_het = [Symbolics.hessian(t, [a, b]) for t in test_exp]
function h_dense_arraymat_het_julia!(out, x)
    a, b, c = x
    out[1] .= [exp(a[1]) * exp(b[1]) exp(a[1]) * exp(b[1]); exp(a[1]) * exp(b[1]) exp(a[1]) * exp(b[1])]
    out[2] .= [0 0; 0 0]
end

h_dense_arraymat_het_str = Symbolics.build_function(h_dense_arraymat_het, [a, b, c])
h_dense_arraymat_het_cse_str = Symbolics.build_function(h_dense_arraymat_het, [a, b, c], cse=true)
h_dense_arraymat_het_ip! = eval(h_dense_arraymat_het_str[2])
h_dense_arraymat_het_cse_ip! = eval(h_dense_arraymat_het_cse_str[2])
out_1_arraymat_het = [Array{Float64}(undef, 2, 2) for i in 1:2]
out_2_arraymat_het = [similar(x) for x in out_1_arraymat_het]
out_3_arraymat_het = [similar(x) for x in out_1_arraymat_het]
h_dense_arraymat_het_julia!(out_1_arraymat_het, input)
h_dense_arraymat_het_ip!(out_2_arraymat_het, input)
h_dense_arraymat_het_cse_ip!(out_3_arraymat_het, input)
@test out_1_arraymat_het == out_2_arraymat_het == out_3_arraymat_het

# Arrays of 1D Vectors
h_dense_arrayvec = [[a, 0, c], [0, 0, 0], [1, a, b]] # same for empty vectors, etc.
function h_dense_arrayvec_julia!(out, x)
    a, b, c = x
    out[1] .= [a[1], 0, c[1]]
    out[2] .= [0, 0, 0]
    out[3] .= [1, a[1], b[1]]
end

h_dense_arrayvec_str = Symbolics.build_function(h_dense_arrayvec, [a, b, c])
h_dense_arrayvec_ip! = eval(h_dense_arrayvec_str[2])
out_1_arrayvec = [Vector{Int64}(undef, 3) for i in 1:3]
out_2_arrayvec = [Vector{Int64}(undef, 3) for i in 1:3]
h_dense_arrayvec_julia!(out_1_arrayvec, input)
h_dense_arrayvec_ip!(out_2_arrayvec, input)
@test out_1_arrayvec == out_2_arrayvec

# Arrays of Arrays of Matrices
h_dense_arrayNestedMat = [[[a 1; b 0], [0 0; 0 0]], [[b 1; a 0], [b c; 0 1]]]
function h_dense_arrayNestedMat_julia!(out, x)
    a, b, c = x
    out[1][1] .= [a[1] 1; b[1] 0]
    out[1][2] .= [0 0; 0 0]
    out[2][1] .= [b[1] 1; a[1] 0]
    out[2][2] .= [b[1] c[1]; 0 1]
end

h_dense_arrayNestedMat_str = Symbolics.build_function(h_dense_arrayNestedMat, [a, b, c])
h_dense_arrayNestedMat_ip! = eval(h_dense_arrayNestedMat_str[2])
out_1_arrayNestedMat = [[rand(Int64, 2, 2), rand(Int64, 2, 2)], [rand(Int64, 2, 2), rand(Int64, 2, 2)]] # avoid undef broadcasting issue
out_2_arrayNestedMat = [[rand(Int64, 2, 2), rand(Int64, 2, 2)], [rand(Int64, 2, 2), rand(Int64, 2, 2)]]
h_dense_arrayNestedMat_julia!(out_1_arrayNestedMat, input)
h_dense_arrayNestedMat_ip!(out_2_arrayNestedMat, input)
@test out_1_arrayNestedMat == out_2_arrayNestedMat

# Arrays of Arrays of Matrices, Heterogeneous element types
test_exp = [exp(a) * exp(b), a]
h_dense_arrayNestedMat_het = [[Symbolics.hessian(t, [a, b]) for t in test_exp], [Num[0 0; 0 0], Num[0 0; 0 0]]]
function h_dense_arrayNestedMat_het_julia!(out, x)
    a, b, c = x
    out[1][1] .= [exp(a[1]) * exp(b[1]) exp(a[1]) * exp(b[1]); exp(a[1]) * exp(b[1]) exp(a[1]) * exp(b[1])]
    out[1][2] .= [0 0; 0 0]
    out[2][1] .= [0 0; 0 0]
    out[2][2] .= [0 0; 0 0]
end
h_dense_arrayNestedMat_het_str = Symbolics.build_function(h_dense_arrayNestedMat_het, [a, b, c])
h_dense_arrayNestedMat_het_ip! = eval(h_dense_arrayNestedMat_het_str[2])
out_1_arrayNestedMat_het = [[rand(Int64, 2, 2), rand(Int64, 2, 2)], [rand(Int64, 2, 2), rand(Int64, 2, 2)]] # avoid undef broadcasting issue
out_2_arrayNestedMat_het = [[rand(Int64, 2, 2), rand(Int64, 2, 2)], [rand(Int64, 2, 2), rand(Int64, 2, 2)]]
h_dense_arrayNestedMat_julia!(out_1_arrayNestedMat_het, input)
h_dense_arrayNestedMat_ip!(out_2_arrayNestedMat_het, input)
@test out_1_arrayNestedMat_het == out_2_arrayNestedMat_het

# ===== Sparse tests =====
# Array of Matrices
h_sparse_arraymat = sparse.([[a 1; b 0], [0 0; 0 0], [a c; 1 0]])
function h_sparse_arraymat_julia!(out, x)
    a, b, c = x
    out[1][1, 1] = a[1]
    out[1][1, 2] = 1
    out[1][2, 1] = b[1]
    out[2] = sparse([0 0; 0 0]) # no undef constructor for SparseMatrixCSC
    out[3][1, 1] = a[1]
    out[3][1, 2] = c[1]
    out[3][2, 1] = 1
end

h_sparse_arraymat_str = Symbolics.build_function(h_sparse_arraymat, [a, b, c])
h_sparse_arraymat_ip! = eval(h_sparse_arraymat_str[2])
h_sparse_arraymat_sparsity_patterns = map(get_sparsity_pattern, h_sparse_arraymat)
out_1_arraymat = [similar(h) for h in h_sparse_arraymat_sparsity_patterns]
out_2_arraymat = [similar(h) for h in h_sparse_arraymat_sparsity_patterns] # can't do similar() because it will just be #undef, with the wrong sparsity pattern
h_sparse_arraymat_ip!(out_2_arraymat, input)
h_sparse_arraymat_julia!(out_1_arraymat, input)
@test out_1_arraymat == out_2_arraymat

# Array of 1D Vectors
h_sparse_arrayvec = sparse.([[a, 0, c], [0, 0, 0], [1, a, b]])
function h_sparse_arrayvec_julia!(out, x)
    a, b, c = x
    out[1][1] = a[1]
    out[1][3] = c[1]
    out[2] = sparse([0, 0, 0]) # necessary because sparsity pattern is 3 elements with 0 stored, not 0 elements
    out[3][1] = 1
    out[3][2] = a[1]
    out[3][3] = b[1]
end

h_sparse_arrayvec_str = Symbolics.build_function(h_sparse_arrayvec, [a, b, c])
h_sparse_arrayvec_ip! = eval(h_sparse_arrayvec_str[2])
h_sparse_arrayvec_sparsity_patterns = map(get_sparsity_pattern, h_sparse_arrayvec)
out_1_arrayvec = [similar(h) for h in h_sparse_arrayvec_sparsity_patterns]
out_2_arrayvec = [similar(h) for h in h_sparse_arrayvec_sparsity_patterns]
h_sparse_arrayvec_julia!(out_1_arrayvec, input)
h_sparse_arrayvec_ip!(out_2_arrayvec, input)
@test out_1_arrayvec == out_2_arrayvec

# Arrays of Arrays of Matrices
h_sparse_arrayNestedMat = [sparse.([[a 1; b 0], [0 0; 0 0]]), sparse.([[b 1; a 0], [b c; 0 1]])]
function h_sparse_arrayNestedMat_julia!(out, x)
    a, b, c = x
    out[1][1][1, 1] = a[1]
    out[1][1][1, 2] = 1
    out[1][1][2, 1] = b[1]
    out[1][2] = sparse([0 0; 0 0])
    out[2][1][1, 1] = b[1]
    out[2][1][1, 2] = 1
    out[2][1][2, 1] = a[1]
    out[2][2][1, 1] = b[1]
    out[2][2][1, 2] = c[1]
    out[2][2][2, 2] = 1
end

h_sparse_arrayNestedMat_str = Symbolics.build_function(h_sparse_arrayNestedMat, [a, b, c])
h_sparse_arrayNestedMat_ip! = eval(h_sparse_arrayNestedMat_str[2])
h_sparse_arrayNestedMat_sparsity_patterns = [map(get_sparsity_pattern, h) for h in h_sparse_arrayNestedMat]
out_1_arrayNestedMat = [[similar(h_sub) for h_sub in h] for h in h_sparse_arrayNestedMat_sparsity_patterns]
out_2_arrayNestedMat = [[similar(h_sub) for h_sub in h] for h in h_sparse_arrayNestedMat_sparsity_patterns]
h_sparse_arrayNestedMat_julia!(out_1_arrayNestedMat, input)
h_sparse_arrayNestedMat_ip!(out_2_arrayNestedMat, input)
@test out_1_arrayNestedMat == out_2_arrayNestedMat

# Additional Tests
# Returning 0-element structures (corresponding to empty Jacobians)
# Arrays of Matrices
h_empty = [[a b; c 0], Array{Num,2}(undef, 0,0)]
h_empty_str = Symbolics.build_function(h_empty, [a, b, c])
h_empty_ip! = eval(h_empty_str[2])
out = [Matrix{Int64}(undef, 2, 2), Matrix{Int64}(undef, 0, 0)]
h_empty_ip!(out, input) # should just not fail

# Array of Vectors
h_empty_vec = [[a, b, c, 0], Vector{Num}(undef,0)]
h_empty_vec_str = Symbolics.build_function(h_empty_vec, [a, b, c])
h_empty_vec_ip! = eval(h_empty_vec_str[2])
out = [Vector{Int64}(undef, 4), Vector{Int64}(undef, 0)]
h_empty_vec_ip!(out, input) # should just not fail

# Arrays of Arrays of Matrices
h_emptyNested = [[[a b; c 0]], Array{Array{Num, 2}}(undef, 0)] # emptyNested array of arrays
h_emptyNested_str = Symbolics.build_function(h_emptyNested, [a, b, c])
h_emptyNested_ip! = eval(h_emptyNested_str[2])
out = [[[1 2;3 4]], Array{Array{Int64,2},1}(undef, 0)]
h_emptyNested_ip!(out, input) # should just not fail
