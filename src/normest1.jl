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
    I = eye(n)
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
        ind = integer(h_ind[:,2])#round(Int,h_ind[:,2])
        if maximum(h) <= dot(Z[:,ind_best],X[:,ind_best])
            itcount = k
            break
        end
        if debug
            println(h[1])
        end
        for j in 1:t
            X[:,j] = I[:,ind[j]]
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
  est_old = 0 # initial estimate
  ind = zeros(n,1) 
  S = zeros(Int8,n,t)  
  S_old = zeros(n,t)  
  int_hist = [] #integer vector recording indices of used unit vectors
  for k=1:itmax+1
    Y = A*X # TODO: make Y preallocated, and in place multiplication
    est, est_indx = findmax([norm(slice(X,1:n,col),1) for col in 1:t])
    if est > est_old || k==2
      ind_best = est_indx
      w = Y[:,ind_best]
    end
    if k>=2 && est <= est_old
      est = est_old
      break
    end
    est_old = est
    copy!(S_old,S)
    if k>itmax
        break
    end
    copy!(S,sign(Y))
    # TODO: If every column of S is parallel to a column of S_old, break
    if t>1
      # TODO:
      # Ensure that no column of S is parallel to another column of S
      # or to a column of Sold by replacing columns of S by rand{−1, 1}.
    end
    Z = transpose(A)*S
    h = Float64[norm(Z[i,:],Inf) for i in 1:n]
    if k>=2 and maximum(h) == h[ind_best]
        break
    end
    # TODO: Sort h in decreasing order, and reorder ind correspondingly
    if t>1
        if ind[1:t] in ind_hist
            break
        end
        # TODO: replace ind[1:t] with the first t indices in
        #       ind[1:n] that are not in ind_hist
    end
    for j in 1:t
        X[:,j] = 
    end
    ind_hist = [int_hist; ind[1:t]]
  end
  return est,
end

function test(k;d=10)
    A = randn(d,d); 
    return ht21(A), ht22(A,k), norm(A,1)
#     res = zeros(k,5)
#     for i in 1:k
#         A = randn(d,d)
#         res[i,1] = norm(A,1)
#         est21, it21 = ht21(A)
#         res[i,2] = est21
#         res[i,3] = it21
#         est22, it22 = ht22(A,d)
#         res[i,4] = est22
#         res[i,5] = it22
#     end
#     res
end
