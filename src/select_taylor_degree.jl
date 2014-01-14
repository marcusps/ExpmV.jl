#function  [M,mv,alpha,unA] = ...
#           select_taylor_degree(A,b,m_max,p_max,prec,shift,bal,force_estm)

function select_taylor_degree(A,
                              b;
                              m_max = 55,
                              p_max = 8,
                              prec = "double",
                              shift,
                              # bal,
                              force_estm = false)


    #SELECT_TAYLOR_DEGREE   Select degree of Taylor approximation.
    #   [M,MV,alpha,unA] = SELECT_TAYLOR_DEGREE(A,m_max,p_max) forms a matrix M
    #   for use in determining the truncated Taylor series degree in EXPMV
    #   and EXPMV_TSPAN, based on parameters m_max and p_max.
    #   MV is the number of matrix-vector products with A or A^* computed.
    
    #   Reference: A. H. Al-Mohy and N. J. Higham, Computing the action of
    #   the matrix exponential, with an application to exponential
    #   integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.
    
    #   Awad H. Al-Mohy and Nicholas J. Higham, October 26, 2010.
    
    if p_max < 2 || m_max > 60 || m_max + 1 < p_max*(p_max - 1)
        error('>>> Invalid p_max or m_max.')
    end
    
    n = size(A,1);
    
    # if bal
    #     [D B] = balance(A);
    #     if norm(B,1) < norm(A,1), A = B; end
    # end
    
    if prec == "double" 
        load theta_taylor
    elseif prec == "single"
        load theta_taylor_single
    elseif prec == "half"
        load theta_taylor_half
    end
    
    if shift
        mu = trace(A)/n;
        A = A-mu*speye(n);
    end
    
    mv = 0;
    
    if !force_estm
        normA = norm(A,1)
    end
    
    if !force_estm && normA <= 4*theta(m_max)*p_max*(p_max + 3)/(m_max*size(b,2));
        unA = 1;
        c = normA;
        alpha = c*ones(p_max-1,1);
    else
        unA = 0;
        eta = zeros(p_max,1); 
        alpha = zeros(p_max-1,1);
        for p = 1:p_max
            (c,k) = normAm(A,p+1);
            c = c^(1/(p+1));
            mv = mv + k;
            eta[p] = c;
        end
        for p = 1:p_max-1
            alpha[p] = max(eta(p),eta(p+1));
        end
    end
    
    M = zeros(m_max,p_max-1);
    
    for p = 2:p_max
        for m = p*(p-1)-1 : m_max
            M[m,p-1] = alpha(p-1)/theta(m);
        end
    end
    
    return (M,mv,alpha,unA)
    
end