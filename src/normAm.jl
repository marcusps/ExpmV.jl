function normAm(A,m)
    #NORMAM   Estimate of 1-norm of power of matrix.
    #   NORMAM(A,m) estimates norm(A^m,1).
    #   If A has nonnegative elements the estimate is exact.
    #   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
    #   matrix-vector products computed involving A or A^*.
    
    #   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    #   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
    #   970-989, 2009.
    
    #   Awad H. Al-Mohy and Nicholas J. Higham, September 7, 2010.
    
    t = 1; # Number of columns used by NORMEST1.
    
    n = length(A);

    if isequal(A,abs(A))
        e = ones(n,1);
        for j=1:m         # for positive matrices only
            e = A'*e;
        end
        c = norm(e,Inf);
        mv = m;
    else
        (c,v,w,it) = normest1(@afun_power,t);
        mv = it(2)*t*m;
    end
    
    return (c,mv)
    
end

function afun_power(flag,X)
    #AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.
    
    if flag == "dim"
        Z = n;
    elseif flag == "real"
        if isreal(A)
            Z = 1
        else
            Z = 0
        end
    else
        (p,q) = size(X);
        if p != n, 
            error("Dimension mismatch") 
        end
        
        if flag=="notransp"
            for i = 1:m
                X = A*X; 
            end
        elseif flag=="transp"
            for i = 1:m
                X = A'*X
            end
        end
        
        Z = X;
        
    end
    
    return Z
end
