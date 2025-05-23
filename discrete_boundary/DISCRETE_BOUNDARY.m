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

discrete_boundary_value_fun = @discrete_boundary_value_fvalue;
discrete_boundary_value_grad = @discrete_boundary_value_grad;
discrete_boundary_value_hess = @discrete_boundary_value_hess;

discrete_boundary_value_grad_fd = @discrete_boundary_value_grad_fd;
discrete_boundary_value_hess_fd = @discrete_boundary_value_hess_fd;
    
% ======================= MODIFIED NEWTON ===========================

tic;
for p=1:length(d)

    fprintf('Sto stampando risultati per p = %d\n', p);

    n = 10^d(p);

    % Backtracking parameters
    rho = 0.5;
    c = 10e-4;

    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-8;
    btmax = 20;

    % con il mio x_bar
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

    [xk2_prec, fk2_prec, gradfk_norm2_prec, k2_prec, xseq2_prec, btseq2_prec] = ...
        modified_newton_bcktrck_preconditioning(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
        discrete_boundary_value_grad , discrete_boundary_value_hess, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton_discrete_boundary_value_prec = xk2_prec;

    % [xk2_fd, fk2_fd, gradfk_norm2_fd, k2_fd, xseq2_fd, btseq2_fd] = ...
    %     modified_newton_bcktrck(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
    %     discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_discrete_boundary_value_fd = xk2_fd;

    % [xk2_fd_prec, fk2_fd_prec, gradfk_norm2_fd_prec, k2_fd_prec, xseq2_fd_prec, btseq2_fd_prec] = ...
    %     modified_newton_bcktrck_preconditioning(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
    %     discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
    %     kmax, tolgrad, c1, rho, btmax);
    % x_newton_discrete_boundary_value_fd_prec = xk2_fd_prec;

    fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk2, k2, gradfk_norm2);
    fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk2_prec, k2_prec, gradfk_norm2_prec);
    % fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk2_fd, k2_fd, gradfk_norm2_fd);
    % fprintf("n = %d | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, fk2_fd_prec, k2_fd_prec, gradfk_norm2_fd_prec);

    % con i 10 punti generati uniformemente in un ipercubo
    for i = 1:num_points

        x0_i = x_bar_discrete_boundary_value + 2 * rand(n,1) - 1;

        [xk2, fk2, gradfk_norm2, k2, xseq2, btseq2] = ...
            modified_newton_bcktrck(x0_i, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad , discrete_boundary_value_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        x_newton_discrete_boundary_value = xk2;
    
        [xk2_prec, fk2_prec, gradfk_norm2_prec, k2_prec, xseq2_prec, btseq2_prec] = ...
            modified_newton_bcktrck_preconditioning(x0_i, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad , discrete_boundary_value_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        x_newton_discrete_boundary_value_prec = xk2_prec;

        % [xk2_fd, fk2_fd, gradfk_norm2_fd, k2_fd, xseq2_fd, btseq2_fd] = ...
        %     modified_newton_bcktrck(x0_i, discrete_boundary_value_fun, ...
        %     discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
        %     kmax, tolgrad, c1, rho, btmax);
        % x_newton_discrete_boundary_value_fd = xk2_fd;
        % 
        % [xk2_fd_prec, fk2_fd_prec, gradfk_norm2_fd_prec, k2_fd_prec, xseq2_fd_prec, btseq2_fd_prec] = ...
        %     modified_newton_bcktrck_preconditioning(x0_i, discrete_boundary_value_fun, ...
        %     discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
        %     kmax, tolgrad, c1, rho, btmax);
        % x_newton_discrete_boundary_value_fd_prec = xk2_fd_prec;

        fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
            n, i, fk2, k2, gradfk_norm2);
        fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
            n, i, fk2_prec, k2_prec, gradfk_norm2_prec);
        % fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
        %     n, i, fk2_fd, k2_fd, gradfk_norm2_fd);
        % fprintf("n = %d | Punto #%d | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
        %     n, i, fk2_fd_prec, k2_fd_prec, gradfk_norm2_fd_prec);

    end
end
tempoNewton = toc; 
disp(['Tempo trascorso: ', num2str(tempoNewton), ' secondi']);

% ======================= NELDER-MEAD ===========================

tic;
for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    fprintf('Sto stampando risultati per n = %d\n', n);

    % con il mio x_bar
    x_bar_discrete_boundary_value = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar_discrete_boundary_value(i) = i*h*(1-i*h);
    end

    simplex_discrete_boundary_value = nelder_mead_n(x_bar_discrete_boundary_value, discrete_boundary_value_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

    % Restituisco valore migliore di ogni simplesso
    simplex_discrete_boundary_value = simplex_discrete_boundary_value(:,1);
    
    fprintf("Nelder-Mead | n = %d | #%d | f(x) = %.4e\n", n, i, discrete_boundary_value_fun(simplex_discrete_boundary_value));

    for i = 1:num_points

        x0_i = x_bar_discrete_boundary_value + 2 * rand(n,1) - 1;
    
        simplex_i = nelder_mead_n(x0_i, discrete_boundary_value_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    
        x_best_i = simplex_i(:,1);
    
        fprintf("Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, discrete_boundary_value_fun(x_best_i));

    end
    
end
tempoNM = toc;
disp(['Tempo trascorso: ', num2str(tempoNM), ' secondi']);


