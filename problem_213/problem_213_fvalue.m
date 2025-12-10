% Computes the function value for the Problem 213

function F = problem_213_fvalue(x)

    n = length(x);

    %f = zeros(n,1);
    h = 1/(n+1);

 
    x = [0; x; 1];

    xkm1 = x(1:n);
    xk = x(2:n+1);
    xkp1 = x(3:n+2);

    % for k=1:n
    % 
    %     xkm1 = x(k);
    %     xk = x(k+1);
    %     xkp1 = x(k+2);
    % 
    %     f(k)= 2*xk + (h^2)*(xk + sin(xk))-xkm1-xkp1;
    % 
    % end
    f = 2*xk + h^2*(xk + sin(xk)) - xkm1 - xkp1;

    F = 0.5*sum(f.^2);

end