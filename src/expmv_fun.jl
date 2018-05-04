export expmv

function expmv(t::Number, A, b; M = [], prec = "double", shift = false, full_term = false, prnt = false)

    n = size(A, 1)

    if shift
        mu = trace(A)/n
        #mu = full(mu) # Much slower without the full!
        A = A-mu*speye(n)
    end

    if isempty(M)
        tt = 1
        (M,alpha,unA) = select_taylor_degree(t*A,b)

    else
        tt = t

        mvd = 0
    end

    tol =
      if prec == "double"
          2.0^(-53)
      elseif prec == "single"
          2.0^(-24)
      elseif prec == "half"
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

    # if prnt
    #     fprintf("m = %2.0f, s = %g, m_actual = ", m, s)
    # end

    for i = 1:s
        c1 = norm(b,Inf);
        for k = 1:m
            b = (t/(s*k))*(A*b);
            f =  f + b;
            c2 = norm(b,Inf);
            if !full_term
                if c1 + c2 <= tol*norm(f,Inf)
                    # if prnt
                    #     fprintf(" %2.0f, ", k)
                    # end
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
