function F = problem_213_fvalue(x)

    n = length(x);

    f = zeros(n,1);
    h = 1/(n+1);

    x_first = 0;
    x_last = 1;
    x = [x_first; x; x_last];

    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
            
        f(k)= 2*xk + (h^2)*(xk + sin(xk))-xkm1-xkp1;

    end

    F = 0.5*sum(f.^2);

end