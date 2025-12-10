% Computes the gradient of the function for the Problem 213

function grad = problem_213_grad(x)

    n = length(x);
    h = 1/(n+1);
    
    x = [0; x; 1];

    xkm1 = x(1:n);
    xk   = x(2:n+1);
    xkp1 = x(3:n+2);

    fk = 2*xk + (h^2)*(xk + sin(xk))-xkm1-xkp1;
    df_dxk   = 2 + h^2 * (1+cos(xk));
    
    grad = zeros(n,1);
   % contributo da fk * df/dx_k
    grad = grad + fk .* df_dxk;

    % contributo da fk * df/dx_{k+1} = fk * (-1)
    grad(1:n-1) = grad(1:n-1) - fk(2:n);

    % contributo da fk * df/dx_{k-1} = fk * (-1)
    grad(2:n) = grad(2:n) - fk(1:n-1);

    % for k=1:n
    % 
    %     xkm1 = x(k);
    %     xk = x(k+1);
    %     xkp1 = x(k+2);
    % 
    %     fk = 2*xk + (h^2)*(xk + sin(xk))-xkm1-xkp1;
    % 
    %     df_dxkm1 = -1;
    %     df_dxk   = 2 + h^2 * (1+cos(xk));
    %     df_dxkp1 = -1;
    % 
    %     if k > 1
    %         grad(k-1) = grad(k-1) + (fk * df_dxkm1);
    %     end
    %     grad(k) = grad(k) + (fk * df_dxk);
    %     if k < n
    %         grad(k+1) = grad(k+1) + (fk * df_dxkp1);
    %     end
    % 
    % end

end
