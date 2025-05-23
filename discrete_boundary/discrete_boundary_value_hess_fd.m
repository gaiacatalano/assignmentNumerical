function H = discrete_boundary_value_hess_fd(x, hstep)
    % Approximate Hessian for the discrete BVP via finite differences
    %
    % Input:
    %   x     - n x 1 vector
    %   hstep - finite difference step (default: 1e-5)
    %
    % Output:
    %   H     - n x n sparse Hessian approximation

    if nargin < 2
        hstep = 1e-5;
    end

    n = length(x);
    H = sparse(n, n);
    h = 1 / (n + 1);

    for i = 1:n
        % Perturb positively
        xp = x;
        xp(i) = xp(i) + hstep;
        gp = discrete_boundary_value_grad_fd(xp, h);

        % Perturb negatively
        xm = x;
        xm(i) = xm(i) - hstep;
        gm = discrete_boundary_value_grad_fd(xm, h);

        % Approximate column j of the Hessian
        H(:,i) = (gp - gm) / (2 * hstep);
    end

    % Enforce symmetry
    H = 0.5 * (H + H');
end
