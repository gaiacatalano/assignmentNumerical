function H = chained_rosenbrock_hess_fd(x, hstep)
    if nargin < 2
        hstep = 1e-5;
    end

    n = length(x);
    H = sparse(n, n);

    for i = 1:n
        % Perturb x in direction j
        xp = x;
        xp(i) = xp(i) + hstep;
        gp = chained_rosenbrock_grad_fd(xp);  % funzione strutturata

        xm = x;
        xm(i) = xm(i) - hstep;
        gm = chained_rosenbrock_grad_fd(xm);

        % Central difference approximation of ∂grad/∂x_j
        H(:, i) = (gp - gm) / (2 * hstep);
    end

    % Enforce symmetry
    H = 0.5 * (H + H');


end