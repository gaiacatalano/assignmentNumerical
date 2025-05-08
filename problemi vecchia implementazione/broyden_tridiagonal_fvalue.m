function F = broyden_tridiagonal_fvalue(x)

    n = length(x);

    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];
    
    f = zeros(n,1);

    for k=1:n
        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);

        f(k) = (3-2*xk)*xk - xkm1 -2*xkp1 + 1;
    end

    F = 0.5 * sum(f.^2);

end

