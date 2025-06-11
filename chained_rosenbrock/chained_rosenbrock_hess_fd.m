function H = chained_rosenbrock_hess_fd(x, h)
    if nargin < 2
        h = 1e-5;
    end

    n = length(x);
    d0 = zeros(n,1);
    dp1 = zeros(n,1);
    dm1 = zeros(n,1);

    for i = 1:n
        xi = x(i);        
        if i<n
            xip1 = x(i+1);
            dm1(i) = -200*(2*xi + h);
            dp1(i+1) = -200*(2*xi + h);
            d0(i) = 400*(xi^2 - xip1) + 2 + 200*(h^2 + 4*xi^2);
        end
        if i ~= 1
            d0(i)= d0(i) + 200;
            
        end    
    end

    H = spdiags([dm1 d0 dp1], [-1 0 1], n, n);

    % Enforce symmetry
    H = 0.5 * (H + H');


end