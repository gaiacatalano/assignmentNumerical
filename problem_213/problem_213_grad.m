% Computes the gradient of the function for the Problem 213

function grad = problem_213_grad(x)

    n = length(x);
    h = 1/(n+1);
    
    x_first = 0;
    x_last = 1;
    x = [x_first; x; x_last];

    grad = zeros(n,1);

    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
            
        fk = 2*xk + (h^2)*(xk + sin(xk))-xkm1-xkp1;

        df_dxkm1 = -1;
        df_dxk   = 2 + h^2 * (1+cos(xk));
        df_dxkp1 = -1;

        if k > 1
            grad(k-1) = grad(k-1) + (fk * df_dxkm1);
        end
        grad(k) = grad(k) + (fk * df_dxk);
        if k < n
            grad(k+1) = grad(k+1) + (fk * df_dxkp1);
        end

    end

end
