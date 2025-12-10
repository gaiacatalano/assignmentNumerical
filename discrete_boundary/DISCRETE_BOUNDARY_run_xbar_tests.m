clc; clear; close all;

addpath('../');
rng(349131);

d = 3;         % n = 10^3 = 1000
n = 10^d;

tolgrad = 1e-6;
kmax    = 10000;
rho     = 0.5;
c1      = 1e-4;
btmax   = 1000;

% Functions
f_fun  = @discrete_boundary_value_fvalue;
f_grad = @discrete_boundary_value_grad;
f_hess = @discrete_boundary_value_hess;

f_grad_fd = @discrete_boundary_value_grad_fd;
f_hess_fd = @discrete_boundary_value_hess_fd;

% Output file
fid = fopen('output_xbar.txt','w');

fprintf('Running tests with starting point x_bar...\n');
fprintf(fid, '===== Starting point: x_bar, n = %d =====\n\n', n);

% Construct x_bar
h = 1/(n+1);
i = (1:n)';
x_bar = i*h .* (1 - i*h);

% Sweep on hstep
for k = 2:2:24
    fprintf('k = %d\n', k);

    if k > 12
        hstep = 10^(-(k-12));
        hstep_i = 1;
    else
        hstep = 10^(-k);
        hstep_i = 0;
    end

    if hstep < 1e-10
        fprintf(fid, "Skipping FD for hstep = %.1e (too small)\n", hstep);
        do_fd = false;
    else
        do_fd = true;
    end

    % --- Newton analitico ---
    tic;
    [xk, fk, gk, kiter] = modified_newton_backtracking(x_bar, f_fun, f_grad, f_hess, ...
                               kmax, tolgrad, c1, rho, btmax);
    t_analytic = toc;

    % --- Newton analitico + precond ---
    tic;
    [xk_prec, fk_prec, gk_prec, kiter_prec] = ...
        modified_newton_backtracking_preconditioning(x_bar, f_fun, f_grad, f_hess, ...
                              kmax, tolgrad, c1, rho, btmax);
    t_analytic_prec = toc;

    if do_fd
        % --- Newton FD ---
        tic;
        [xk_fd, fk_fd, gk_fd, kiter_fd] = modified_newton_backtracking(x_bar, f_fun, ...
                        f_grad_fd, f_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        t_fd = toc;
    
        % --- Newton FD + precond ---
        tic;
        [xk_fd_prec, fk_fd_prec, gk_fd_prec, kiter_fd_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar, f_fun, ...
                        f_grad_fd, f_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        t_fd_prec = toc;
    end

    

    % Write to file
    fprintf(fid, "===== hstep = %.2e  (k = %d) =====\n", hstep, k);
    fprintf(fid, "Newton analytic | t = %.4f | iter = %d | f = %.3e | grad = %.3e\n", ...
        t_analytic, kiter, fk, gk);
    fprintf(fid, "Newton analytic + prec | t = %.4f | iter = %d | f = %.3e | grad = %.3e\n", ...
        t_analytic_prec, kiter_prec, fk_prec, gk_prec);
    if do_fd
        fprintf(fid, "Newton FD | t = %.4f | iter = %d | f = %.3e | grad = %.3e\n", ...
            t_fd, kiter_fd, fk_fd, gk_fd);
        fprintf(fid, "Newton FD + prec | t = %.4f | iter = %d | f = %.3e | grad = %.3e\n\n", ...
            t_fd_prec, kiter_fd_prec, fk_fd_prec, gk_fd_prec);
    end
end

fclose(fid);
fprintf('Saved results to output_xbar.txt\n');
