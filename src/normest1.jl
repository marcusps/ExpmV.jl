using Compat

function ht21(A;itmax=20)
    n,m = size(A)
    if n!=m
        error("norm 1 estimation can only be applied to square matrices")
    end
    x = ones(n)/n
    γ = 0
    itcount = itmax
    for k in 1:itmax
        y = A*x
        ξ = sign(y)
        z = A'*ξ
        if norm(z,Inf)<=dot(z,x) && k>1
            γ=norm(y,1)
            itcount = k
            break
        end
        x = zeros(Float64,n); x[indmax(abs(z))] = 1.
    end
    b = Float64[ (-1)^(i+1) * (1+(i-1)/(n-1)) for i in 1:n]
    return max(γ,norm(A*b,1)/norm(b,1)), itcount
end

function ht22(A,t::Int;itmax=10,debug=false)
    n,m = size(A)
    if n!=m
        error("norm 1 estimation can only be applied to square matrices")
    end
    # start with a random matrix of normalized columns
    X = randn(n,t)
    for col=1:t
        scale!(slice(X,1:n,col),1/norm(slice(X,1:n,col),1))
    end
    Id = eye(n)
    h = Float64[]
    itcount = itmax
    for k in 1:itmax
        Y = A*X
        g = Float64[norm(Y[:,j],1) for j in 1:t]
        ind_best = indmax(g)
        sort!(g,rev=true)
        if debug
            println(g[1])
        end
        S = sign(Y)
        Z = A'*S
        h = Float64[maximum(abs(vec(Z[i,:]))) for i in 1:n]
        ind = [1:n;]
        h_ind = sortrows(hcat(h,ind),by=x->x[1],rev=true)
        h   = h_ind[:,1]
        ind = @compat round(Integer,h_ind[:,2])
        if maximum(h) <= dot(Z[:,ind_best],X[:,ind_best])
            itcount = k
            break
        end
        if debug
            println(h[1])
        end
        for j in 1:t
            X[:,j] = Id[:,ind[j]]
        end
    end
    return h[1], itcount
end

function ht23(A,t::Int;itmax=2)
    n,m = size(A)
    if n!=m
        error("norm 1 estimation can only be applied to square matrices")
    end
    # start with a random matrix of normalized columns
    X = randn(size(A,1),t)
    for col=1:t
        scale!(slice(X,1:n,col),1/norm(slice(X,1:n,col),1))
    end
    Id = eye(n)
    est_old = 0 # initial old estimate
    est = 0 # initial estimate
    ind = zeros(n,1) 
    S = zeros(Int8,n,t)  
    S_old = zeros(n,t)  
    ind_hist = [] # integer vector recording indices of used unit vectors
    itcount = itmax
    for k=1:itmax+1
        Y = A*X # TODO: make Y preallocated, and in place multiplication
        est, est_indx = findmax([norm(slice(Y,1:n,col),1) for col in 1:t])
        if est > est_old || k==2
            ind_best = est_indx
            w = Y[:,ind_best]
        end
        if k>=2 && est <= est_old
            est = est_old
            itcount = k
            break
        end
        est_old = est
        copy!(S_old,S)
        if k>itmax
            break
        end
        copy!(S,sign(Y))
        # If every column of S is parallel to a column of S_old, break
        if parallel_cols(S,S_old)
            itcount = k
            break
        end
        if t>1
            randomize_parallel!(S,S_old)
        end
        Z = transpose(A)*S
        h = Float64[norm(Z[i,:],Inf) for i in 1:n]
        if k>=2 && maximum(h) == h[ind_best]
            itcount = k
            break
        end
        # Sort h in decreasing order, and reorder ind correspondingly
        ind = [1:n;]
        h_ind = sortrows(hcat(h,ind),by=x->x[1],rev=true)
        h   = h_ind[:,1]
        #ind = integer(h_ind[:,2]) 
        ind = @compat round(Integer,h_ind[:,2])
        if t>1
            if ind[1:t] ⊆ ind_hist
                itcount = k
                break
            end
            # replace ind[1:t] with the first t indices in
            # ind[1:n] that are not in ind_hist
            ind[1:t] = filter(x->!(x ⊆ ind_hist),ind)[1:t]
        end
        for j in 1:t
            X[:,j] = Id[:,ind[j]]
        end
        ind_hist = [ind_hist; ind[1:t]]
    end
    return est, itcount
end

function parallel_cols(S,So)
    n = size(S,1)
    # the assumption is that S,So are matrices with elements ±1,
    # so being parallel would translate to abs value of inner product ≈ n
    return any(map(x->isapprox(x,n),abs(So'*S)),1) |> all
end

function randomize_parallel!(S,So)
    n = size(S,1)
    # Ensure no column of S is parallel to a column of S
    ips = map(x->isapprox(x,n),abs(S'*S))
    for i=1:n
        ips[i,i] = false
        if any(ips[:,i])
            S[:,i] = rand([-1,1],n)
        end
    end
    # Ensure no column of S is parallel to a column of So
    ips = map(x->isapprox(x,n),abs(So'*S))
    for i=1:n
        ips[i,i] = false
        if any(ips[:,i])
            S[:,i] = rand([-1,1],n)
        end
    end
end

function test(k;d=10)
    A = randn(d,d); 
    return ht21(A), ht22(A,k), ht23(A,k), norm(A,1)
end
