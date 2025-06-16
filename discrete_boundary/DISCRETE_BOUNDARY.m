clc
clear
close all

% Aggiungo quello che sta nella cartella fuori al path
addpath('../'); 

% Seed
rng(349131);

% Dimension  
d = 3:1:3; 
num_points = 10;

% Stopping parameters
tol = 1e-6;
kmax = 800;

% Chiamo le funzioni

discrete_boundary_value_fun = @discrete_boundary_value_fvalue;
discrete_boundary_value_grad = @discrete_boundary_value_grad;
discrete_boundary_value_hess = @discrete_boundary_value_hess;

discrete_boundary_value_grad_fd = @discrete_boundary_value_grad_fd;
discrete_boundary_value_hess_fd = @discrete_boundary_value_hess_fd;

fid = fopen('output_discrete_boundary.txt', 'w');
    
% ======================= MODIFIED NEWTON ===========================
for p=1:length(d)

    fprintf('Sto stampando risultati per p = %d\n', p);

    n = 10^d(p);
    fprintf(fid, "n = %d\n", n);

    % Backtracking parameters
    rho = 0.5;

    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-4;
    btmax = 800; 

    % con il mio x_bar
    x_bar_discrete_boundary_value = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar_discrete_boundary_value(i) = i*h*(1-i*h);
    end

    hstep_i=0;

    for k = 2:2:12
        if k > 12
            hstep  = 10^(-(k-12));
            hstep_i = 1;
        else
            hstep  = 10^(-k);
        end

        tic;
        [xk2, fk2, gradfk_norm2, k2, xseq2, btseq2] = ...
            modified_newton_bcktrck(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad , discrete_boundary_value_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn = toc;
        x_newton_discrete_boundary_value = xk2;
    
        tic;
        [xk2_prec, fk2_prec, gradfk_norm2_prec, k2_prec, xseq2_prec, btseq2_prec] = ...
            modified_newton_bcktrck_preconditioning(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad , discrete_boundary_value_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn_prec = toc;
    
        x_newton_discrete_boundary_value_prec = xk2_prec;
    
        tic;
        [xk2_fd, fk2_fd, gradfk_norm2_fd, k2_fd, xseq2_fd, btseq2_fd] = ...
            modified_newton_bcktrck(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd = toc;
        x_newton_discrete_boundary_value_fd = xk2_fd;
    
        tic;
        [xk2_fd_prec, fk2_fd_prec, gradfk_norm2_fd_prec, k2_fd_prec, xseq2_fd_prec, btseq2_fd_prec] = ...
            modified_newton_bcktrck_preconditioning(x_bar_discrete_boundary_value, discrete_boundary_value_fun, ...
            discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd_prec = toc;
        x_newton_discrete_boundary_value_fd_prec = xk2_fd_prec;
         
        fprintf(fid, "Tempo di esecuzione Modified Newton: %.4f\n", tempo_mn);
        fprintf(fid, "n = %d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_discrete_boundary_value), fk2, k2, gradfk_norm2);
        fprintf(fid, "Tempo di esecuzione Modified Newton con precondizionamento: %.4f\n", tempo_mn_prec);
        fprintf(fid, "n = %d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_discrete_boundary_value_prec), fk2_prec, k2_prec, gradfk_norm2_prec);
        fprintf(fid, "Tempo di esecuzione Modified Newton con differenze finite: %.4f\n", tempo_mn_fd);
        fprintf(fid, "n = %d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_discrete_boundary_value_fd), fk2_fd, k2_fd, gradfk_norm2_fd);
        fprintf(fid, "Tempo di esecuzione Modified Newton con differenze finite e precondizionamento: %.4f\n", tempo_mn_fd_prec);
        fprintf(fid, "n = %d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_discrete_boundary_value_fd_prec), fk2_fd_prec, k2_fd_prec, gradfk_norm2_fd_prec);
    
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
    
            [xk2_fd, fk2_fd, gradfk_norm2_fd, k2_fd, xseq2_fd, btseq2_fd] = ...
                modified_newton_bcktrck(x0_i, discrete_boundary_value_fun, ...
                discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            x_newton_discrete_boundary_value_fd = xk2_fd;
    
            [xk2_fd_prec, fk2_fd_prec, gradfk_norm2_fd_prec, k2_fd_prec, xseq2_fd_prec, btseq2_fd_prec] = ...
                modified_newton_bcktrck_preconditioning(x0_i, discrete_boundary_value_fun, ...
                discrete_boundary_value_grad_fd , discrete_boundary_value_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            x_newton_discrete_boundary_value_fd_prec = xk2_fd_prec;
    
             fprintf(fid, "n = %d | Punto #%d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                 n, i, norm(x_newton_discrete_boundary_value), fk2, k2, gradfk_norm2);
            fprintf(fid, "n = %d | Punto #%d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(x_newton_discrete_boundary_value_prec), fk2_prec, k2_prec, gradfk_norm2_prec);
            fprintf(fid, "n = %d | Punto #%d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(x_newton_discrete_boundary_value_fd), fk2_fd, k2_fd, gradfk_norm2_fd);
            fprintf(fid, "n = %d | Punto #%d | norma 2 di x = %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(x_newton_discrete_boundary_value_fd_prec), fk2_fd_prec, k2_fd_prec, gradfk_norm2_fd_prec);
    
        end
    end
end


% ======================= NELDER-MEAD ===========================

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

    tic;
    simplex_discrete_boundary_value = nelder_mead_n(x_bar_discrete_boundary_value, discrete_boundary_value_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    tempo_nelder_mead = toc;

    % Restituisco valore migliore di ogni simplesso
    simplex_discrete_boundary_value = simplex_discrete_boundary_value(:,1);
    
    fprintf(fid, "Tempo di esecuzione Nelder Mead: %.4f\n", tempo_nelder_mead);
    fprintf("Nelder-Mead | n = %d | #%d | f(x) = %.4e\n", n, i, discrete_boundary_value_fun(simplex_discrete_boundary_value));

    for i = 1:num_points

        x0_i = x_bar_discrete_boundary_value + 2 * rand(n,1) - 1;
    
        simplex_i = nelder_mead_n(x0_i, discrete_boundary_value_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    
        x_best_i = simplex_i(:,1);
    
        fprintf(fid, "Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, discrete_boundary_value_fun(x_best_i));

    end
    
end

fclose(fid);



