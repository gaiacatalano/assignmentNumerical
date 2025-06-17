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
tol = 1e-5;
kmax = 30000;

% treashold norm_grad
epsilon = 1e-4;
count_success_newton = 0;
count_failure_newton = 0;
count_success_nelder = 0;
count_failure_nelder = 0;

% Chiamo le funzioni

problem_213_fun = @problem_213_fvalue;
problem_213_grad = @problem_213_grad;
problem_213_hess = @problem_213_hess;

problem_213_grad_fd = @problem_213_grad_fd;
problem_213_hess_fd = @problem_213_hess_fd;

fid = fopen('output_problem_213.txt', 'w');

% ======================= MODIFIED NEWTON ===========================

fprintf(fid, "Modified Newton method\n");

for p=1:length(d)

    fprintf('Sto stampando risultati per p = %d\n', p);

    n = 10^d(p);
    %fprintf(fid, "n = %d\n", n);

    % Backtracking parameters
    rho = 0.5;
    c = 1e-4;

    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-8;
    btmax = 60;

    % con il mio x_bar
    x_bar_problem_213 = ones(n,1);

    hstep_i = 0;

    for k = 2:2:24        
             
        if k > 12
            hstep  = 10^(-(k-12));
            hstep_i=1;
        else
            hstep  = 10^(-k);
        end
           
        tic;
        [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
            modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn = toc;
        x_newton_problem_213 = xk3;

        if gradfk_norm3 < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end

        %r_k = log(norm(problem_213_grad(k))/norm(problem_213_grad(k+1))) / log(norm(grad_km1)/norm(grad_km2));

        
        tic;
        [xk3_prec, fk3_prec, gradfk_norm3_prec, k3_prec, xseq3_prec, btseq3_prec] = ...
            modified_newton_bcktrck_preconditioning(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn_prec = toc;
        x_newton_problem_213_prec = xk3_prec;

        if gradfk_norm3_prec < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end
    
        tic;
        [xk3_fd, fk3_fd, gradfk_norm3_fd, k3_fd, xseq3_fd, btseq3_fd] = ...
            modified_newton_bcktrck(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad_fd , problem_213_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd = toc;
        x_newton_problem_213_fd = xk3_fd;

        if gradfk_norm3_fd < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end
    
        tic;
        [xk3_fd_prec, fk3_fd_prec, gradfk_norm3_fd_prec, k3_fd_prec, xseq3_fd_prec, btseq3_fd_prec] = ...
            modified_newton_bcktrck_preconditioning(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad_fd , problem_213_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd_prec = toc;
        x_newton_problem_213_fd_prec = xk3_fd_prec;

        if gradfk_norm3_fd_prec < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end
        
        fprintf(fid, "Tempo di esecuzione Modified Newton: %.4f\n", tempo_mn);
        fprintf(fid, "n = %d |norma2 di x= %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_problem_213), fk3, k3, gradfk_norm3);
        fprintf(fid, "Tempo di esecuzione Modified Newton con precondizionamento: %.4f\n", tempo_mn_prec);
        fprintf(fid, "n = %d |norma2 di x= %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_problem_213_prec), fk3_prec, k3_prec, gradfk_norm3_prec);
        fprintf(fid, "Tempo di esecuzione Modified Newton con differenze finite: %.4f\n", tempo_mn_fd);
        fprintf(fid, "n = %d |norma2 di x= %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_problem_213_fd), fk3_fd, k3_fd, gradfk_norm3_fd);
        fprintf(fid, "Tempo di esecuzione Modified Newton con differenze finite e precondizionamento: %.4f\n", tempo_mn_fd_prec);
        fprintf(fid, "n = %d |norma2 di x= %.2e | f(x) = %.4e | iter = %d | norm grad = %.2e\n", n, norm(x_newton_problem_213_fd_prec), fk3_fd_prec, k3_fd_prec, gradfk_norm3_fd_prec);
    
        % con i 10 punti generati uniformemente in un ipercubo
        for i = 1:num_points
            
            x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;
    
            % Newton classico
            [xk_rand, fk_rand, gradfk_norm_rand, k_rand] = ...
                modified_newton_bcktrck(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);

            if gradfk_norm_rand < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            % Newton precondizionato
            [xk_rand_prec, fk_rand_prec, gradfk_norm_rand_prec, k_rand_prec] = ...
                modified_newton_bcktrck_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);

            if gradfk_norm_rand_prec < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            [xk_rand_fd, fk_rand_fd, gradfk_norm_rand_fd, k_rand_fd, xseq_fd, btseq_fd] = ...
                modified_newton_bcktrck(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            %x_newton_problem_213_fd = xk_fd;

            if gradfk_norm_rand_fd < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            [xk_rand_fd_prec, fk_rand_fd_prec, gradfk_norm_rand_fd_prec, k_rand_fd_prec, xseq_fd_prec, btseq_fd_prec] = ...
                modified_newton_bcktrck_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            %x_newton_problem_213_fd_prec = xk_fd_prec;

            if gradfk_norm_rand_fd_prec < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            fprintf(fid, "n = %d | Punto #%d |norma2 di x= %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(xk_rand), fk_rand, k_rand, gradfk_norm_rand);
            fprintf(fid, "n = %d | Punto #%d |norma2 di x= %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(xk_rand_prec), fk_rand_prec, k_rand_prec, gradfk_norm_rand_prec);
            fprintf(fid, "n = %d | Punto #%d |norma2 di x= %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(xk_rand_fd), fk_rand_fd, k_rand_fd, gradfk_norm_rand_fd);
            fprintf(fid, "n = %d | Punto #%d |norma2 di x= %.2e | f(x) = %.4e | iter = %d (norm grad = %.2e)\n", ...
                n, i, norm(xk_rand_prec), fk_rand_fd_prec, k_rand_fd_prec, gradfk_norm_rand_fd_prec);
    
        end
    end

    fprintf(fid, "n = %d | number of successes of modified newton: %.4f\n", n, count_success_newton);
    fprintf(fid, "n = %d |number of failures of modified newton: %.4f\n",n,  count_failure_newton);

end


% ======================= NELDER-MEAD ===========================

fprintf(fid, "Nelder Mead");

for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    fprintf('Sto stampando i simplessi per n= %d\n', n)

    % con il mio x_bar
    x_bar_problem_213 = ones(n,1);

    tic;
    simplex_problem_213 = nelder_mead_n(x_bar_problem_213, problem_213_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    tempo_nelder_mead = toc; 
    D = max(pdist(simplex_problem_213'));
    if D < epsilon
        count_success_nelder = count_success_nelder +1;
    else
        count_failure_nelder = count_failure_newton + 1;
    end
    

    % Restituisco valore migliore del simplesso
    simplex_problem_213 = simplex_problem_213(:,1);

    fprintf(fid, "Tempo di esecuzione Nelder Mead: %.4f\n", tempo_nelder_mead);
    fprintf(fid, "Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, problem_213_fun(simplex_problem_213));

    for i = 1:num_points
        x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;

        simplex_i = nelder_mead_n(x0_i, problem_213_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

        Di = max(pdist(simplex_i'));
        if Di < epsilon
            count_success_nelder = count_success_nelder +1;
        else
            count_failure_nelder = count_failure_newton + 1;
        end

        x_best_i = simplex_i(:,1);

        fprintf(fid, "Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, problem_213_fun(x_best_i));
    end

    fprintf(fid, "n = %d | number of successes of nelder mead: %.4f\n", n, count_success_nelder);
    fprintf(fid, "n = %d | number of failures of nelder mead: %.4f\n", n, count_failure_nelder);
    
end



fclose(fid);



