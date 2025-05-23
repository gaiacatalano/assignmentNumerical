clc
clear
close all

% Aggiungo quello che sta nella cartella fuori al path
addpath('../'); 

% Seed
rng(349131);

% Dimension
d = 3:1:5; 
num_points = 10;

% Stopping parameters
tol = 1e-5;
kmax = 200;

% Chiamo le funzioni

problem_213_fun = @problem_213_fvalue;
problem_213_grad = @problem_213_grad;
problem_213_hess = @problem_213_hess;

problem_213_grad_fd = @problem_213_grad_fd;
problem_213_hess_fd = @problem_213_hess_fd;


% ======================= MODIFIED NEWTON ===========================

for p=1:length(d)

    disp('Sto stampando risultati per p=', p)

    n = 10^d(p);

    % Backtracking parameters
    rho = 0.5;
    c = 10e-4;

    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-8;
    btmax = 20;

    % con il mio x_bar
    x_bar_problem_213 = ones(n,1);

    [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
        modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
        problem_213_grad , problem_213_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_problem_213 = xk3;

    [xk3_prec, fk3_prec, gradfk_norm3_prec, k3_prec, xseq3_prec, btseq3_prec] = ...
        modified_newton_bcktrck_preconditioning(x_bar_problem_213, problem_213_fun, ...
        problem_213_grad , problem_213_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_problem_213_prec = xk3_prec;

    % [xk3_fd, fk3_fd, gradfk_norm3_fd, k3_fd, xseq3_fd, btseq3_fd] = ...
    %     modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
    %     problem_213_grad_fd , problem_213_hess_fd, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_problem_213_fd = xk3_fd;

    fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk3, k3, gradfk_norm3);
    fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk3_prec, k3_prec, gradfk_norm3_prec);
    
    % con i 10 punti generati uniformemente in un ipercubo
    for i = 1:num_points
        x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;

        % Newton classico
        [xk_rand, fk_rand, gradfk_norm_rand, k_rand] = ...
            modified_newton_bcktrck(x0_i, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);

        % Newton precondizionato
        [xk_rand_prec, fk_rand_prec, gradfk_norm_rand_prec, k_rand_prec] = ...
            modified_newton_bcktrck_preconditioning(x0_i, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);

        fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
            n, i, fk_rand, k_rand, gradfk_norm_rand);
        fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
            n, i, fk_rand_prec, k_rand_prec, gradfk_norm_rand_prec);

    end

end


% ======================= NELDER-MEAD ===========================

for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    fprintf('Sto stampando i simplessi per n= %d\n', n)

    % con il mio x_bar
    x_bar_problem_213 = ones(n,1);

    simplex_problem_213 = nelder_mead_n(x_bar_problem_213, problem_213_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

    % Restituisco valore migliore del simplesso
    simplex_problem_213 = simplex_problem_213(:,1);

    fprintf("Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, problem_213_fun(simplex_problem_213));

    for i = 1:num_points
        x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;

        simplex_i = nelder_mead_n(x0_i, problem_213_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

        x_best_i = simplex_i(:,1);

        fprintf("Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, problem_213_fun(x_best_i));
    end
    
end



