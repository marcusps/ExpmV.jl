using LinearAlgebra
using SparseArrays

function degree_selector(t, M, U, p)
    C = ceil.(abs.(t)*M)'*U
    C = zero_to_inf.(C)

    # idx is a CarthesianIndex if C' is a Matrix, or a scalar if C' is a row
    # vector. idx[1] extract the first index, i.e. row, of the CarthesianIndex
    cost, idx = findmin(C')
    m = idx[1]

    if cost == Inf
        cost = 0
    end
    s = max(cost/m,1);
    return (m, s)
end

function zero_to_inf(x::Number)
    if x == 0.
        Inf
    else
        x
    end
end
