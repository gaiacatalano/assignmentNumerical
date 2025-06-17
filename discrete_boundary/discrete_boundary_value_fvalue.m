% Computes the function value for the Discrete Boundary Value problem

function F = discrete_boundary_value_fvalue(x)

    n = length(x);
    f = zeros(n,1);
    h = 1/(n+1);

    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];

    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
            
        f(k)= (2*xk-xkm1-xkp1 + (h^2 * (xk + k*h + 1)^3)/2)^2;

    end

    F = sum(f.^2);

end

