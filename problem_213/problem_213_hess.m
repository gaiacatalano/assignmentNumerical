% Computes the Hessian of the Problem 213

function hess = problem_213_hess(x)
    
    n = length(x);
    h = 1/(n+1);

    x_first = 0;
    x_last = 1;
    x = [x_first; x; x_last];   

    d0 = zeros(n,1);
    dp1 = zeros(n,1);
    dp2 = zeros(n,1);
    dm1 = zeros(n,1);
    dm2 = zeros(n,1);

    for i = 1:n
        
        fk = 2*x(i+1) + (h^2)*(x(i+1) + sin(x(i+1)))-x(i)-x(i+2);
        d0(i) = (2 + (1 + cos(x(i+1))*h^2))^2 - sin(x(i+1))*h^2*fk;
    
        if i ~= 1
            d0(i) = d0(i) + 1;
        end
            
        if i ~= n
            d0(i) = d0(i) + 1;
            dm1(i) = - 4 - (2 + cos(x(i+1)) + cos(x(i+2)))*h^2;
        end
    
        if i < n-1
            dm2(i) = 1;
        end
        
    end

    hess = spdiags([dm2 dm1 d0 [0; dp1(1:end-1)] [0; 0; dp2(1:end-2)]], [-2 -1 0 1 2], n, n);
    hess = 0.5 * (hess + hess');
    
end