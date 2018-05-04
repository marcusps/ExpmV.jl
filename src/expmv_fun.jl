export expmv

"""
    expmv(t, A, b; <keyword arguments>)

Returns the matrix-vector product ``exp(t A) b`` where `A` is a ``n × n`` sparse
real or complex matrix, `b` is a vector of ``n`` real or complex elements and `t` is
a parameter (or a `StepRangeLen` object representing a range of values).

# Arguments
* `t`: `Number` or `StepRangeLen` object
* `A`: a `n × n` real or complex sparse matrix
* `b`: an `n`-vector
* `M = []`: manually set the degree of the Taylor expansion
* `precision = "double"`: can be `"double"`, `"single"` or `"half"`.
* `shift = false`: set to `true` to apply a shift in order to reduce the norm of A
        (see Sec. 3.1 of the paper)
* `full_term = false`: set to `true` to evaluate the full Taylor expansion instead
        of truncating when reaching the required precision
"""
function expmv(t::Number, A::SparseMatrixCSC, b::Vector; M = nothing,
                precision = "double", shift = false, full_term = false)
    n = size(A, 1)

    if shift
        mu = trace(A)/n
        A = A-mu*speye(n)
    end

    if M == nothing
        tt = 1
        (M,alpha,unA) = select_taylor_degree(t*A,b)
    else
        tt = t
    end

    tol =
      if precision == "double"
          2.0^(-53)
      elseif precision == "single"
          2.0^(-24)
      elseif precision == "half"
          2.0^(-10)
      end

    s = 1;

    if t == 0
        m = 0;
    else
        (m_max,p) = size(M);
        U = diagm(1:m_max);
        C = ( (ceil.(abs.(tt)*M))'*U );
        zero_els = find(x->x==0, C)
        for el in zero_els
            C[el] = Inf
        end
        if p > 1
            cost,m = findmin(minimum(C,1)); # cost is the overall cost.
        else
            cost,m = findmin(C);  # when C is one column. Happens if p_max = 2.
        end
        if cost == Inf
            cost = 0
        end
        s = max(cost/m,1);
    end

    eta = 1;

    if shift
        eta = exp(t*mu/s)
    end

    f = b;

    for i = 1:s
        c1 = norm(b,Inf);
        for k = 1:m
            b = (t/(s*k))*(A*b);
            f =  f + b;
            c2 = norm(b,Inf);
            if !full_term
                if c1 + c2 <= tol*norm(f,Inf)
                    break
                end
                c1 = c2;
            end

        end
        f = eta*f
        b = f
    end

    return f
end
