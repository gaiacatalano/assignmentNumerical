function grad = chained_rosenbrock_grad(x)

    n = length(x);

    grad = zeros(n,1);

    for i = 2:n
        xim1 = x(i-1);
        xi = x(i);

        df_dxim1 = 400*(xim1^2 - xi)*xim1 + 2*(xim1 - 1);
        df_dxi   = -200*(xim1^2 - xi);

        grad(i-1) = grad(i-1) + df_dxim1;
        grad(i) = grad(i) + df_dxi;
    end

end