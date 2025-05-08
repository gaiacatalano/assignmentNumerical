function[F, gradf, Hess] = discrete_boundary_value_problem(n, x_bar)

    grad = zeros(n,1);
    f = zeros(n,1);
    h = 1/(n+1);

    d0 = zeros(n,1);
    dp1 = zeros(n,1);
    dp2 = zeros(n,1);

    x_first = 0;
    x_last = 0;
    % x_bar qui è già un vettore
    x = [x_first; x_bar; x_last];

    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
            
        f(k)= (2*xk-xkm1-xkp1 + (h^2 * (xk + k*h + 1)^3)/2)^2;

        % Gradient
        df_dxkm1 = -1;
        df_dxk   = 2 + 3/2 * h^2 * (xk + k*h + 1)^2;
        df_dxkp1 = -1;

        if k > 1
            grad(k-1) = grad(k-1) + 2*(f(k) * df_dxkm1);
        end
        grad(k) = grad(k) + 2*(f(k) * df_dxk);
        if k < n
            grad(k+1) = grad(k+1) + 2*(f(k) * df_dxkp1);
        end
        
        %Hessian
        % grad_fk = zeros(n,1);
        % if k > 1
        %     grad_fk(k-1) = df_dxkm1; 
        % end
        % grad_fk(k) = df_dxk;
        % if k < n 
        %     grad_fk(k+1) = df_dxkp1; 
        % end
        % 
        % H = H + 2*(grad_fk * grad_fk');
        % H(k,k) = H(k,k) + 2 * (f(k) * (3*h^2 * (xk + k*h + 1)));

        if k > 1
            d0(k) = d0(k) + df_dxkp1^2;
        end
        d0(k) = d0(k) + 2*df_dxk^2 + 2*f(k)*(3*h^2 * (xk + k*h + 1));
        if k < n
            d0(k) = d0(k) + df_dxkm1^2;
        end

        dp1(k) = dp1(k) + 2*(df_dxk*df_dxkp1 + df_dxkm1*df_dxk);
        dp2(k) = dp2(k) + 2*df_dxkp1^2;
    end

    H = spdiags([dp2 dp1 d0 dp1 dp2], [-2 -1 0 1 2], n, n);
    H = 0.5 * (H + H');

    F = @(x) sum(f.^2);
    gradf = @(x) grad;
    Hess = @(x) H;

end