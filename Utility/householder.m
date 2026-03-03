function Atilde = householder(A)
    Atilde = A;
    n = size(A, 2)-1;
    m = size(A,1) - n;
    

    for k = 1:n
        
        sigma = sign(Atilde(k,k) + (Atilde(k,k)==0))*sqrt(sum(Atilde(k:m+n,k).^2));
        if sigma == 0
            continue
        end
        u_i = Atilde(k:m+n,k);
        u_i(1) = u_i(1) + sigma;
        Atilde(k,k) = -sigma;

        beta = 1/(sigma*u_i(1));
        subA = Atilde(k:m+n, k+1:n+1);
        Atilde(k:m+n, k+1:n+1) = subA - (beta * u_i) * (u_i' * subA);
        
        Atilde(k+1:m+n,k) = 0;

    end
end