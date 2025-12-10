% Computes the function value for the Discrete Boundary Value problem

function F = discrete_boundary_value_fvalue(x)

    n = length(x);
    %f = zeros(n,1);
    h = 1/(n+1);

    %x_first = 0;
    %x_last = 0;
    %x = [x_first; x; x_last];
    x_ext = [0; x; 0];

    k = (1:n)';

    % Extract shifted versions of x
    xkm1 = x_ext(1:n);     % x(k-1)
    xk   = x_ext(2:n+1);   % x(k)
    xkp1 = x_ext(3:n+2);   % x(k+1)

    % Compute f_k vectorized
    fk = 2*xk - xkm1 - xkp1 + 0.5*h^2 * (xk + h*k + 1).^3;

    % Compute final objective
    F = (fk.' * fk);   % = sum(fk.^2)

    % for k=1:n
    % 
    %     xkm1 = x(k);
    %     xk = x(k+1);
    %     xkp1 = x(k+2);
    % 
    %     f(k)= (2*xk-xkm1-xkp1 + (h^2 * (xk + k*h + 1)^3)/2)^2;
    % 
    % end
    % 
    % F = sum(f.^2);

end

