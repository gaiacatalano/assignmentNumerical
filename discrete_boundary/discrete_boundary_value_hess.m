% Computes the Hessian of the Discrete Boundary Value problem 

function hess = discrete_boundary_value_hess(x)
    
    n = length(x);
    h = 1/(n+1);

    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];    

    d0 = zeros(n,1);
    dm1 = zeros(n,1);
    dm2 = zeros(n,1);

    for k=1:n

        fk = 2*x(k+1)-x(k)-x(k+2) + (h^2 * (x(k+1) + k*h + 1)^3)/2;
        d0(k) = 2*(2 + 3/2 * h^2 * (x(k+1) + k*h + 1)^2)^2 + 2*fk*(3*h^2 * (x(k+1) + k*h + 1));

        if k ~= 1
            d0(k) = d0(k) + 1;            
        end
        
        if k ~= n
            d0(k) = d0(k) + 1;
            dm1(k) = -(2 + 3/2 * h^2 * (x(k+1) + k*h + 1)^2) - ...
                (2 + 3/2 * h^2 * (x(k+2) + k*h + 1)^2);
        end

        if k < n-1
            dm2(k) = 2;
        end

    end

    hess = spdiags([dm2 dm1 d0 [0; dm1(1:end-1)] [0; 0; dm2(1:end-2)]], [-2 -1 0 1 2], n, n);
    hess = 0.5 * (hess + hess');

end