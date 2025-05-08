function grad = broyden_tridiagonal_grad(x)

    n = length(x);
    
    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];

    grad = zeros(n,1);

    for k = 1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);

        df_dxkm1 = -1;
        df_dxk   = 3 - 4*xk;
        df_dxkp1 = -2;

        fk = (3-2*xk)*xk - xkm1 -2*xkp1 + 1;

        if k > 1
            grad(k-1) = grad(k-1) + fk * df_dxkm1;
        end
        grad(k) = grad(k) + fk * df_dxk;
        if k < n
            grad(k+1) = grad(k+1) + fk * df_dxkp1;
        end

    end
        
end

