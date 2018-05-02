function degree_selector(t, M, U, p)
    C = ceil.(abs.(t)*M)'*U
    C = zero_to_inf.(C)
    if p > 1
        cost, m = findmin(minimum(C,1)); # cost is the overall cost.
    else
        cost, m = findmin(C);  # when C is one column. Happens if p_max = 2.
    end
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
