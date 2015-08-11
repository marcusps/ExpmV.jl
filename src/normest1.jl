function ht22(A,t::Int)
    # TODO: if A is not nxn (square), error
    n,m = size(A)
    if n!=m
        error("norm 1 estimation can only be applied to square matrices")
    end
    # start with a random matrix of normalized columns
    X = randn(size(A,1),t)
    for col=1:t
        scale!(slice(X,1:n,col),1/norm(slice(X,1:n,col),1))
    end
    I = eye(n)
    while true
        Y = A*X
        g = Float64[norm(Y[:,j],1) for j in 1:t]
        _, ind_best = findmax(g)
        sort!(g)
        S = sign(Y)
        Z = A'*S
        h = Float64[norm(Z[i,:],Inf) for i in 1:n]
        ind = [1:n;]
        if maximum(h) <= Z[:,ind_best]'*X[:,ind_best]
            break
        end
        h_ind = sortrows(hcat(h,ind),by=x->x[1])
        h   = h_ind[:,1]
        ind = h_ind[:,2]
        for j in 1:t
            X[:,j] = I[:,ind[j]]
        end
    end
    return h
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

function test(k,d=10)
    res = zeros(k,2)
    for i in 1:k
        A = randn(d,d)+1im*randn(d,d)
        res[i,1] = norm(A,1)
        res[i,2] = ht22(A,6)
    end
    res
end