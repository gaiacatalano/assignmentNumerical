clc; clear; close all;

addpath('../');
rng(349131);

d = 3;
n = 10^d;

tolgrad = 1e-6;
kmax = 10000;
rho  = 0.5;
c1   = 1e-4;
btmax = 1000;
num_points = 10;

% Functions
f_fun  = @discrete_boundary_value_fvalue;
f_grad = @discrete_boundary_value_grad;
f_hess = @discrete_boundary_value_hess;

fid = fopen('output_random.txt', 'w');

fprintf('Running random starting point tests...\n');
fprintf(fid, "===== Random points, n = %d =====\n\n", n);

% Construct x_bar
h = 1/(n+1);
i = (1:n)';
x_bar = i*h .* (1 - i*h);

for k = 1:num_points
    fprintf('Point #%d...\n', k);

    x0 = x_bar + 2*rand(n,1) - 1;

    % Newton analytic
    [xk, fk, gk, kiter] = modified_newton_backtracking(x0, f_fun, f_grad, f_hess, ...
                                kmax, tolgrad, c1, rho, btmax);

    % Newton analytic + precond
    [xk_prec, fk_prec, gk_prec, kiter_prec] = ...
        modified_newton_backtracking_preconditioning(x0, f_fun, f_grad, f_hess, ...
                                kmax, tolgrad, c1, rho, btmax);

    % Write to file
    fprintf(fid, "Point %d | Newton: iter = %d | f = %.3e | grad = %.3e\n", ...
                k, kiter, fk, gk);
    fprintf(fid, "Point %d | Newton + prec: iter = %d | f = %.3e | grad = %.3e\n\n", ...
                k, kiter_prec, fk_prec, gk_prec);
end

fclose(fid);
fprintf('Saved results to output_random.txt\n');
