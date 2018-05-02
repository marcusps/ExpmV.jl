function expmv(t::StepRangeLen, A, b; M = [], prec = "double", shift = true, full_term = false, prnt = false)
              
    t0 = Float64(t.ref)
    tmax = Float64(t.ref + (t.len - 1.) * t.step)
    q = t.len - 1 # For some reason
    n = size(A, 1)

    mu = 0.
    if shift
        mu = trace(A)/n
        #mu = full(mu) # Much slower without the full!
        A = A-mu*speye(n)
    end

    force_estm = false; temp = (tmax -t0)*norm(A, 1);
    if (prec == "single" || prec == "half" && temp > 85.496) || ( prec == "double" && temp > 63.152)
       force_estm = true;
    end

    if isempty(M)
        (M, mvd, alpha, unA) = select_taylor_degree(A, b, prec=prec, shift=shift, force_estm=force_estm)
        mv = mvd
    else
        mv = 0
    end

    tol =
      if prec == "double"
          2.0^(-53)
      elseif prec == "single"
          2.0^(-24)
      elseif prec == "half"
          2.0^(-10)
      end

    X = zeros(n, q+1);
    (m_max, p) = size(M);
    U = diagm(1:m_max);

    temp, s = degree_selector(tmax - t0, M, U, p)
    h = (tmax - t0)/q;

    X[:,1], mvnew = expmvi(t0,A,b,M=M,prec=prec,shift=shift,prnt=prnt);
    mv += mvnew
    d = Int(floor(q/s))
    j = Int(floor(q/d))
    r = q - d * j
    z = X[:,1]

    m_opt, = degree_selector(d, M, U, p)
    dr = d
    for i = 1:j+1
        if i > j
            dr = r
        end
        K = zeros(n, m_opt+1)
        K[:,1] = z
        m = 0
        for k = 1:dr
            f = z
            c1 = norm(z, Inf)
            for p = 1:m_opt
                if p > m
                    K[:,p+1] = (h/p)*(A*K[:,p])
                    mv += 1
                end

                temp = (k^p)*K[:,p+1]
                f += temp
                c2 = norm(temp, Inf)

                if c1 + c2 <= tol * norm(f, Inf)
                    break
                end
                c1 = c2
            end

            m = max(m,p)

            X[:, k + (i - 1) * d + 1] = exp(k * h * mu) * f
        end

        if i <= j
            z = X[:,i*d+1]
        end
    end

    return X
end