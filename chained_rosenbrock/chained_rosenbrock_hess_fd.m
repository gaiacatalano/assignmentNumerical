function H = chained_rosenbrock_hess_fd(x, h, bool_hstep_i)
    if nargin < 2
        h = 1e-5;
    end

    if bool_hstep_i==1
        hstep_i = abs(x)*h;
    end

    n = length(x);
    d0 = zeros(n,1);
    dp1 = zeros(n,1);
    dm1 = zeros(n,1);

    for i = 1:n

        if bool_hstep_i==1
            h = hstep_i(i);
        end

        %xi = x(i);        
        if i ~= n
            %xip1 = x(i+1);
            dm1(i) = -200*(2*x(i) + h);
            %dp1(i+1) = -200*(2*xi + h);
            d0(i) = 400*(x(i)^2 - x(i+1)) + 2 + 200*(h^2 + 4*x(i)^2);
        end
        if i ~= 1
            d0(i)= d0(i) + 200;
            
        end    
    end

    H = spdiags([dm1, d0, [0; dp1(1:end-1)]], [-1 0 1], n, n);

    % Enforce symmetry
    H = 0.5 * (H + H');


end