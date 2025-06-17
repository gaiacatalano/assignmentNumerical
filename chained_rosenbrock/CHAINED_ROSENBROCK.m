clc
clear
close all

addpath('../'); 

% Seed
rng(349131);

d = 3:1:5; 
num_points = 10;

% Stopping parameters
tol = 1e-06;
kmax = 10000;

% Function recall
chained_rosenbrock_fun = @chained_rosenbrock_fvalue;
chained_rosenbrock_grad = @chained_rosenbrock_grad;
chained_rosenbrock_hess = @chained_rosenbrock_hess;

chained_rosenbrock_grad_fd = @chained_rosenbrock_grad_fd;
chained_rosenbrock_hess_fd = @chained_rosenbrock_hess_fd;

fid = fopen('output_chained_rosenbrock.txt', 'w');
    
% ======================= MODIFIED NEWTON ===========================

fprintf(fid, 'Modified Newton Method\n');
for p=1:length(d)

    fprintf('Displaying results for p = %d\n', p);

    n = 10^d(p);

    % Newton parameters
    tolgrad = 1e-06;
    rho = 0.5;
    c1 = 1e-04;
    btmax = 1000;

    % First starting point x_bar
    x_bar_chained_rosenbrock = zeros(n,1);         
    for i = 1:n
        if mod(i, 2) == 1
            x_bar_chained_rosenbrock(i) = -1.2;
        else
            x_bar_chained_rosenbrock(i) = 1.0;
        end
    end 

    hstep_i=0;

    for k = 2:2:24
        if k > 12
            hstep  = 10^(-(k-12));
            hstep_i = 1;
        else
            hstep = 10^(-k);
        end

        tic;
        [xk, fk, gradfk_norm, k, xseq, btseq] = ...
            modified_newton_backtracking(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
            chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
        tempo_mn = toc;
        x_newton_chained_rosenbrock = xk';
    
        tic;
        [xk_prec, fk_prec, gradfk_norm_prec, k_prec, xseq_prec, btseq_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
            chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
        tempo_mn_prec = toc;
        x_newton_chained_rosenbrock_prec = xk';
    
        tic;
        [xk_fd, fk_fd, gradfk_norm_fd, k_fd, xseq_fd, btseq_fd] = ...
            modified_newton_backtracking(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad_fd , ...
            chained_rosenbrock_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd = toc;
        x_newton_chained_rosenbrock_fd = xk_fd;
    
        tic;
        [xk_fd_prec, fk_fd_prec, gradfk_norm_fd_prec, k_fd_prec, xseq_fd_prec, btseq_fd_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar_chained_rosenbrock, chained_rosenbrock_fun, chained_rosenbrock_grad_fd , ...
            chained_rosenbrock_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd_prec = toc;
        x_newton_chained_rosenbrock_fd_prec = xk_fd_prec;
    
        fprintf(fid, "Modified Newton execution time: %.4f\n", tempo_mn);
        fprintf(fid, "n = %d | f(x) = %.4e | iter = %d | grad norm = %.2e\n", n, fk, k, gradfk_norm);
        fprintf(fid, "Modified Newton execution time with preconditioning: %.4f\n", tempo_mn_prec);
        fprintf(fid, "f(x) = %.4e | iter = %d | grad norm = %.2e\n", fk_prec, k_prec, gradfk_norm_prec);
        fprintf(fid, "Modified Newton execution time with finite differences: %.4f\n", tempo_mn_fd);
        fprintf(fid, "f(x) = %.4e | iter = %d | grad norm = %.2e\n", fk_fd, k_fd, gradfk_norm_fd);
        fprintf(fid, "Modified Newton execution time with finite differences and preconditioning: %.4f\n", tempo_mn_fd_prec);
        fprintf(fid, "f(x) = %.4e | iter = %d | grad norm = %.2e\n", fk_fd_prec, k_fd_prec, gradfk_norm_fd_prec);
        
        % With 10 starting points generated with uniform distribution in a hyper-cube
        for i = 1:num_points
    
            x0_i = x_bar_chained_rosenbrock + 2 * rand(n,1) - 1;
    
            [xk, fk, gradfk_norm, k, xseq, btseq] = ...
                modified_newton_backtracking(x0_i, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
                chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
            x_newton_chained_rosenbrock = xk;
    
            [xk_prec, fk_prec, gradfk_norm_prec, k_prec, xseq_prec, btseq_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, chained_rosenbrock_fun, chained_rosenbrock_grad , ...
                chained_rosenbrock_hess, kmax, tolgrad, c1, rho, btmax);
            x_newton_chained_rosenbrock_prec = xk;
    
            [xk_fd, fk_fd, gradfk_norm_fd, k_fd, xseq_fd, btseq_fd] = ...
                modified_newton_backtracking(x0_i, chained_rosenbrock_fun, chained_rosenbrock_grad_fd , ...
                chained_rosenbrock_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            x_newton_chained_rosenbrock_fd = xk_fd;
    
            [xk_fd_prec, fk_fd_prec, gradfk_norm_fd_prec, k_fd_prec, xseq_fd_prec, btseq_fd_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, chained_rosenbrock_fun, chained_rosenbrock_grad_fd , ...
                chained_rosenbrock_hess_fd, kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            x_newton_chained_rosenbrock_fd_prec = xk_fd_prec;
    
            fprintf(fid, "n = %d | Point #%d | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, fk, k, gradfk_norm);
            fprintf(fid, "n = %d | Point #%d | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, fk_prec, k_prec, gradfk_norm_prec);
            fprintf(fid, "n = %d | Point #%d | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, fk_fd, k_fd, gradfk_norm_fd);
            fprintf(fid, "n = %d | Point #%d | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, fk_fd_prec, k_fd_prec, gradfk_norm_fd_prec);
    
        end

    end

end


% ======================= NELDER-MEAD ===========================

fprintf(fid, "Nelder Method\n");

for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    fprintf('Displaying results for n = %d\n', n);

    % First starting point x_bar
    x_bar_chained_rosenbrock = zeros(n,1);        
    for i = 1:n
        if mod(i, 2) == 1
            x_bar_chained_rosenbrock(i) = -1.2;
        else
            x_bar_chained_rosenbrock(i) = 1.0;
        end
    end 

    tic;
    simplex_chained_rosenbrock = nelder_mead_n(x_bar_chained_rosenbrock, chained_rosenbrock_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    tempo_nelder_mead = toc;

    simplex_chianed_rosenbrock = simplex_chained_rosenbrock(:,1);

    fprintf(fid, "Nelder-Mead execution time: %.4f\n", tempo_nelder_mead);
    fprintf(fid, "Nelder-Mead | n = %d | #%d | f(x) = %.4e\n", n, i, chained_rosenbrock_fun(simplex_chianed_rosenbrock));
    
    % With 10 starting points generated with uniform distribution in a hyper-cube
    for i = 1:num_points

        x0_i = x_bar_chained_rosenbrock + 2 * rand(n,1) - 1;
        
        simplex_i = nelder_mead_n(x0_i, chained_rosenbrock_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    
        x_best_i = simplex_i(:,1);

        fprintf(fid, "Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, chained_rosenbrock_fun(x_best_i));

     end
    
end

fclose(fid);


