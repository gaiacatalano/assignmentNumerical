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

% Treashold norm_grad
epsilon = 1e-4;
count_success_newton = 0;
count_failure_newton = 0;
count_success_nelder = 0;
count_failure_nelder = 0;

% Function recall
problem_213_fun = @problem_213_fvalue;
problem_213_grad = @problem_213_grad;
problem_213_hess = @problem_213_hess;

problem_213_grad_fd = @problem_213_grad_fd;
problem_213_hess_fd = @problem_213_hess_fd;

fid = fopen('output_problem_213.txt', 'w');

% ======================= MODIFIED NEWTON ===========================

fprintf(fid, "Modified Newton method\n");

for p=1:length(d)

    fprintf('Displaying results for p = %d\n', p);

    n = 10^d(p);
    fprintf(fid, "n = %d\n", n);

    % Newton parameters
    tolgrad = 1e-06;
    rho = 0.5;
    c1 = 1e-4;
    btmax = 1000;

    % First starting point x_bar
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
            modified_newton_backtracking(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn = toc;
        x_newton_problem_213 = xk3;

        if gradfk_norm3 < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end
        
        tic;
        [xk3_prec, fk3_prec, gradfk_norm3_prec, k3_prec, xseq3_prec, btseq3_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar_problem_213, problem_213_fun, ...
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
            modified_newton_backtracking(x_bar_problem_213, problem_213_fun, ...
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
            modified_newton_backtracking_preconditioning(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad_fd , problem_213_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd_prec = toc;
        x_newton_problem_213_fd_prec = xk3_fd_prec;

        if gradfk_norm3_fd_prec < epsilon
            count_success_newton = count_success_newton + 1;
        else
            count_failure_newton = count_failure_newton + 1;
        end
        
        fprintf(fid, "Execution time Modified Newton: %.4f\n", tempo_mn);
        fprintf(fid, "n = %d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d | grad norm = %.2e\n", n, norm(x_newton_problem_213), fk3, k3, gradfk_norm3);
        fprintf(fid, "Execution time Modified Newton with preconditioning: %.4f\n", tempo_mn_prec);
        fprintf(fid, "n = %d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d | grad norm = %.2e\n", n, norm(x_newton_problem_213_prec), fk3_prec, k3_prec, gradfk_norm3_prec);
        fprintf(fid, "Execution time Modified Newton with finite differences: %.4f\n", tempo_mn_fd);
        fprintf(fid, "n = %d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d | grad norm = %.2e\n", n, norm(x_newton_problem_213_fd), fk3_fd, k3_fd, gradfk_norm3_fd);
        fprintf(fid, "Execution time Modified Newton with finite differences and preconditioning: %.4f\n", tempo_mn_fd_prec);
        fprintf(fid, "n = %d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d | grad norm = %.2e\n", n, norm(x_newton_problem_213_fd_prec), fk3_fd_prec, k3_fd_prec, gradfk_norm3_fd_prec);

        % With 10 starting points generated with uniform distribution in a hyper-cube
        for i = 1:num_points
            
            x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;
    
            % Newton 
            [xk_rand, fk_rand, gradfk_norm_rand, k_rand] = ...
                modified_newton_backtracking(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);

            if gradfk_norm_rand < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            % Newton with preconditioning
            [xk_rand_prec, fk_rand_prec, gradfk_norm_rand_prec, k_rand_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);

            if gradfk_norm_rand_prec < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            [xk_rand_fd, fk_rand_fd, gradfk_norm_rand_fd, k_rand_fd, xseq_fd, btseq_fd] = ...
                modified_newton_backtracking(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);

            if gradfk_norm_rand_fd < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            [xk_rand_fd_prec, fk_rand_fd_prec, gradfk_norm_rand_fd_prec, k_rand_fd_prec, xseq_fd_prec, btseq_fd_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);

            if gradfk_norm_rand_fd_prec < epsilon
                count_success_newton = count_success_newton + 1;
            else
                count_failure_newton = count_failure_newton + 1;
            end
    
            fprintf(fid, "n = %d | Point #%d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, norm(xk_rand), fk_rand, k_rand, gradfk_norm_rand);
            fprintf(fid, "n = %d | Point #%d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, norm(xk_rand_prec), fk_rand_prec, k_rand_prec, gradfk_norm_rand_prec);
            fprintf(fid, "n = %d | Point #%d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
                n, i, norm(xk_rand_fd), fk_rand_fd, k_rand_fd, gradfk_norm_rand_fd);
            fprintf(fid, "n = %d | Point #%d | 2-norm of x = %.2e | f(x) = %.4e | iter = %d (grad norm = %.2e)\n", ...
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

    fprintf('Displaying results for n = %d\n', n);

    % First starting point x_bar
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
    
    simplex_problem_213 = simplex_problem_213(:,1);

    fprintf(fid, "Nelder-Mead execution time: %.4f\n", tempo_nelder_mead);
    fprintf(fid, "Nelder-Mead | n=%d | #%d | f(x)=%.4e\n", n, i, problem_213_fun(simplex_problem_213));

    % With 10 starting points generated with uniform distribution in a hyper-cube
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



