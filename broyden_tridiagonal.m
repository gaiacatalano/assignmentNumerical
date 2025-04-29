function[F, gradf, Hess] = broyden_tridiagonal(n, x_bar)

    grad = zeros(n,1);
    H = zeros(n,n);
    f = zeros(n,1);

    x_first = 0;
    x_last = 0;
    x = x_bar*ones(n,1);
    x = [x_first; x; x_last];
   
    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);

        f(k) = (3-2*xk)*xk - xkm1 -2*xkp1 + 1;
        
        % Gradient
        df_dxkm1 = -1;
        df_dxk   = 3 - 4*xk;
        df_dxkp1 = -2;

        if k > 1
            grad(k-1) = grad(k-1) + f(k) * df_dxkm1;
        end
        grad(k) = grad(k) + f(k) * df_dxk;
        if k < n
            grad(k+1) = grad(k+1) + f(k) * df_dxkp1;
        end
        
        %Hessian
        grad_fk = zeros(n,1);
        if k > 1
            grad_fk(k-1) = df_dxkm1; 
        end
        grad_fk(k) = df_dxk;
        if k < n 
            grad_fk(k+1) = df_dxkp1; 
        end
        
        H = H + grad_fk * grad_fk';
        H(k,k) = H(k,k) + f(k) * (-4); 
        

    end

    F = @(x) 0.5 * sum(f.^2);
    gradf = @(x) grad;
    Hess = @(x) H;

end