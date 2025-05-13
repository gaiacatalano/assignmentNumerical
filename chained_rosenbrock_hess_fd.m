function H = chained_rosenbrock_hess_fd(x, h)
    if nargin < 2
        h = 1e-5;
    end

    n = length(x);
    %H = sparse(n, n);
    d0 = zeros(n,1);
    d1 = zeros(n,1);

    for i = 2:n
        xi = x(i);
        xim1 = x(i-1);
        % % Perturb x in direction j
        % xp = x;
        % xp(i) = xp(i) + hstep;
        % gp = chained_rosenbrock_grad_fd(xp);  % funzione strutturata
        % 
        % xm = x;
        % xm(i) = xm(i) - hstep;
        % gm = chained_rosenbrock_grad_fd(xm);
        % 
        % % Central difference approximation of ∂grad/∂x_j
        % H(:, i) = (gp - gm) / (2 * hstep);

        if i < n
            xip1 = x(i+1);
            d0(i) = 200*(1+(xi^2 - xip1)*2) + 2 + 200*(h^2 + 4*xi^2);
            %j = i+1;
            d1(i) = -200*(2*xi + h);
        else
            d0(i) = 200; % 100 + 200*(xim1^2 - xi)/h;
        end
    
    end

    H = spdiags([d1 d0 d1], [-1 0 -1], n, n);

    % Enforce symmetry
    H = 0.5 * (H + H');


end