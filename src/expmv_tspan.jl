using LinearAlgebra
using SparseArrays

function expmv(t::StepRangeLen, A, b::AbstractVecOrMat;
                M = nothing, precision = "double", shift = false)

    t0 = Float64(t.ref)
    tmax = Float64(t.ref + (t.len - 1.) * t.step)
    q = t.len - 1
    n = size(A, 1)

    if shift == true && !hasmethod(tr, typeof(A))
        shift = false
        @warn "Shift set to false as $(typeof(A)) doesn't support tr"
    end

    force_estm = !hasmethod(opnorm, Tuple{typeof(A), typeof(1)})
    temp = (tmax - t0) * normAm(A, 1);

    if (precision == "single" || precision == "half" && temp > 85.496) || ( precision == "double" && temp > 63.152)
       force_estm = true;
    end

    if M == nothing
        (M, alpha, unA) = select_taylor_degree(A, b, precision=precision, shift=shift, force_estm=force_estm)
    end

    tol =
      if precision == "double"
          2.0^(-53)
      elseif precision == "single"
          2.0^(-24)
      elseif precision == "half"
          2.0^(-10)
      end

    X = zeros(eltype(A), n, q+1);
    (m_max, p) = size(M);
    U = diagm(0 => 1:m_max);

    temp, s = degree_selector(tmax - t0, M, U, p)
    h = (tmax - t0)/q;

    X[:,1] = expmv(t0,A,b,M=M,precision=precision,shift=shift);

    mu = 0.
    if shift
        mu = tr(A)/n
        A = A - mu * I
    end

    d = max(1, Integer(floor(q/s)))
    j = Integer(floor(q/d))
    r = q - d * j
    z = X[:,1]
    m_opt, = degree_selector(d, M, U, p)
    dr = d
    K = zeros(eltype(A), n, m_opt+1)

    temp = similar(b)
    f = similar(b)

    for i = 1:j+1
        if i > j
            dr = r
        end
        K[:,1] = z
        m = 0
        for k = 1:dr
            f = copy(z)
            c1 = norm(z, Inf)

            p = 1
            for outer p = 1:m_opt
                if p > m
                    @views K[:,p+1] .= (h/p).*(A*K[:,p])
                end

                @views temp .= (Float64(k)^p).*K[:,p+1]
                f .+= temp
                c2 = norm(temp, Inf)
                if c1 + c2 <= tol * norm(f, Inf)
                    break
                end
                c1 = c2
            end

            m = max(m,p)
            @views X[:, k + (i - 1) * d + 1] .= exp(k * h * mu) * f
        end

        if i <= j
            z = X[:,i*d+1]
        end
    end
    return X
end
