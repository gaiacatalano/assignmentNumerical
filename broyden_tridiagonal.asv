function[F] = broyden_tridiagonal(n, x_bar)
    
    x_first = 0;
    x_last = 0;
    x = x_bar*ones(n,1);
    x = [x_first; x; x_last];

    f = zeros(n,1);
    f(1) = (3-2*x(1))*x(1) - x_first -2*x(2) + 1;
    f(n) = (3-2*x(n))*x(n) - x(n-1) -2*x_last + 1;

    for k=2:n-1
       f(k) = (3-2*x(k))*x(k) - x(k-1) -2*x(k+1) + 1;
    end

    F = @(x) 0.5*sum(f.^2);

end