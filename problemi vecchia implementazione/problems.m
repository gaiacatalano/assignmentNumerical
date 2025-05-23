clc
clear
close all

% Seed
rng(349131);

% Dimension
d = 3:1:5; 
num_points = 10;

% Stopping parameters
tol = 1e-5;
kmax = 200;

% Chiamo le funzioni

chained_rosenbrock_fun = @chained_rosenbrock_fvalue;
chained_rosenbrock_grad = @chained_rosenbrock_grad;
chained_rosenbrock_hess = @chained_rosenbrock_hess;

chained_rosenbrock_grad_fd = @chained_rosenbrock_grad_fd;
chained_rosenbrock_hess_fd = @chained_rosenbrock_hess_fd;

discrete_boundary_value_fun = @discrete_boundary_value_fvalue;
discrete_boundary_value_grad = @discrete_boundary_value_grad;
discrete_boundary_value_hess = @discrete_boundary_value_hess;

discrete_boundary_value_grad_fd = @discrete_boundary_value_grad_fd;
discrete_boundary_value_hess_fd = @discrete_boundary_value_hess_fd;

% broyden_tridiagonal_fun = @broyden_tridiagonal_fvalue;
% broyden_tridiagonal_grad = @broyden_tridiagonal_grad;
% broyden_tridiagonal_hess = @broyden_tridiagonal_hess;

problem_213_fun = @problem_213_fvalue;
problem_213_grad = @problem_213_grad;
problem_213_hess = @problem_213_hess;

problem_213_grad_fd = @problem_213_grad_fd;
problem_213_hess_fd = @problem_213_hess_fd;

% Ciclo for per valori di n adatti a Newton
for p=1:length(d)

    disp('Sto stampando risultati per')
    p

    n = 10^d(p);

    % Backtracking parameters
    rho = 0.5;
    c = 10e-4;

    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-8;
    btmax = 20;

    % % Starting points generati random
    % starting_points = 2*rand(n, num_points) - 1;  % genera valori in [-1, 1]
    % starting_points = starting_points + x_bar;    % trasla il centro in x

    % Problem 1
    x_bar_chained_rosenbrock = zeros(n,1);         
    for i = 1:n
        if mod(i, 2) == 1
            x_bar_chained_rosenbrock(i) = -1.2;
        else
            x_bar_chained_rosenbrock(i) = 1.0;
        end
    end 

    [xk, fk, gradfk_norm, k, xseq, btseq] = ...
        modified_newton_bcktrck(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
        chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
    x_newton_chained_rosenbrock = xk;


    [xk, fk, gradfk_norm, k, xseq, btseq] = ...
        modified_newton_bcktrck_preconditioning(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
        chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
    x_newton_chained_rosenbrock_prec = xk;

    % [xk_fd, fk_fd, gradfk_norm_fd, k_fd, xseq_fd, btseq_fd] = ...
    %     modified_newton_bcktrck(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad_fd , ...
    %     chained_rosenbrock_hess_fd, kmax, tolgrad, c1, rho, btmax);
    % x_newton_chained_rosenbrock_fd = xk_fd;

    disp('Sto stampando risultati per Chained')


    % Problem 14
    x_bar_discrete_boundary_value = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar_discrete_boundary_value(i) = i*h*(1-i*h);
    end

    [xk2, fk2, gradfk_norm2, k2, xseq2, btseq2] = ...
        modified_newton_bcktrck(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
        discrete_boundary_value_grad , discrete_boundary_value_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_discrete_boundary_value = xk2;

    [xk2, fk2, gradfk_norm2, k2, xseq2, btseq2] = ...
        modified_newton_bcktrck_preconditioning(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
        discrete_boundary_value_grad , discrete_boundary_value_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_discrete_boundary_value_prec = xk2;

    % [xk2_fd, fk2_fd, gradfk_norm2_fd, k2_fd, xseq2_fd, btseq2_fd] = ...
    %     modified_newton_bcktrck(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
    %     discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_discrete_boundary_value_fd = xk2_fd;
    
    disp('Sto stampando risultati per Discrete')


    % % Problem 31
    % x_broyden_tridiagonal = -1*ones(n,1);
    % 
    % [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    %     modified_newton_bcktrck(x_broyden_tridiagonal, broyden_tridiagonal_fun, broyden_tridiagonal_grad , broyden_tridiagonal_hess, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_broyden_tridiagonal = xk;
    % disp('Sto stampando risultati per Broyden')

    % Problem 83
    x_bar_problem_213 = ones(n,1);

    [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
        modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
        problem_213_grad , problem_213_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_problem_213 = xk3;

    [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
        modified_newton_bcktrck_preconditioning(x_bar_problem_213, problem_213_fun, ...
        problem_213_grad , problem_213_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_problem_213_prec = xk3;

    % [xk3_fd, fk3_fd, gradfk_norm3_fd, k3_fd, xseq3_fd, btseq3_fd] = ...
    %     modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
    %     problem_213_grad_fd , problem_213_hess_fd, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_problem_213_fd = xk3_fd;

    disp('Sto stampando risultati per Problem 213')


end


% Ciclo for per valori di n adatti a Nelder-Mead
for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    % % Starting points generati random
    % starting_points = 2*rand(n, num_points) - 1;  % genera valori in [-1, 1]
    % starting_points = starting_points + x_bar;    % trasla il centro in x
    disp('Sto stampando i simplessi per')
    n

    % Problema 1
    x_bar_chained_rosenbrock = zeros(n,1);        
    for i = 1:n
        if mod(i, 2) == 1
            x_bar_chained_rosenbrock(i) = -1.2;
        else
            x_bar_chained_rosenbrock(i) = 1.0;
        end
    end 

    simplex_chained_rosenbrock = nelder_mead_n(x_bar_chained_rosenbrock, chained_rosenbrock_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    disp('Stampo simplesso chained')


    % Problem 14
    x_bar_discrete_boundary_value = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar_discrete_boundary_value(i) = i*h*(1-i*h);
    end

    simplex_discrete_boundary_value = nelder_mead_n(x_bar_discrete_boundary_value, discrete_boundary_value_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    disp('Stampo simplesso discrete')

    % % Problem 31
    % x_bar_broyden_tridiagonal = -1*ones(n,1);
    % 
    % simplex_broyden_tridiagonal = nelder_mead_n(x_bar_broyden_tridiagonal, broyden_tridiagonal_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    % disp('Stampo simplesso broyden')

    % Problem 83
    x_bar_problem_213 = ones(n,1);

    simplex_problem_213 = nelder_mead_n(x_bar_problem_213, problem_213_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);


    % Restituisco valore migliore di ogni simplesso
    simplex_chianed_rosenbrock = simplex_chained_rosenbrock(:,1);
    simplex_discrete_boundary_value = simplex_discrete_boundary_value(:,1);
    simplex_problem_213 = simplex_problem_213(:,1);
    % simplex_broyden_tridiagonal = simplex_broyden_tridiagonal(:,1);
end



