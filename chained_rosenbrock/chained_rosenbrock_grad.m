function grad = chained_rosenbrock_grad(x)

    n = length(x);

    grad = zeros(n,1);

    for i = 2:n-1
        % xim1 = x(i-1);
        % xi = x(i);
        % 
        % df_dxim1 = 400*(xim1^2 - xi)*xim1 + 2*(xim1 - 1);
        % df_dxi   = -200*(xim1^2 - xi);

        % df_dxim1 = 400*(xi^2 - x(i+1))*xi + 2*(xi - 1);
        % df_dxi   = -200*(xim1^2 - xi);
        % 
        % grad(i-1) = grad(i-1) + df_dxim1;
        % grad(i) = grad(i) + df_dxi;
        grad(i) = 400*x(i)*(x(i)^2 - x(i+1)) + 2*(x(i) - 1) - 200*(x(i-1)^2 - x(i));

      
    end

    grad(1) =  400*x(1)*(x(1)^2 - x(2)) + 2*(x(1) - 1);
    grad(n) =  -200*(x(n-1)^2 - x(n));

end