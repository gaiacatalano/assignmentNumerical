clc
clear
close all

addpath('../'); 

% Seed
rng(349131);

d = 3:1:3; 
num_points = 10;

% Stopping parameters
tol = 1e-06;
kmax = 1000;

% Treashold norm_grad
epsilon = 1e-3;
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
% Newton parameters
tolgrad = 1e-06;
rho = 0.5;
c1 = 1e-4;
btmax = 1000;

fprintf(fid, "Modified Newton method - Tables\n\n");

for p=1:length(d)

    n = 10^d(p);
    fprintf('Running tests for n = %d\n', n);
    fprintf(fid, "============================================================\n");
    fprintf(fid, "n = %d\n", n);
    fprintf(fid, "============================================================\n\n");

    % First starting point x_bar
    x_bar_problem_213 = ones(n,1);

    hstep_i = 0;

    for k = 2 %:2:24        
             
        if k > 12
            hstep  = 10^(-(k-12));
            hstep_i=1;
        else
            hstep  = 10^(-k);
        end
           
        % ============================================================
        % Build a table for THIS (n, hstep)
        % ============================================================
        rows = struct('start',{},'method',{},'iter',{},'xnorm',{},'fval',{},'gnorm',{},'time',{});

        % ------------------ Starting point: x_bar -------------------
        startLabel = 'x_bar';

        % 1) Modified Newton
        tic;
        [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
            modified_newton_backtracking(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn = toc;
        rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton", k3, norm(xk3), fk3, gradfk_norm3, tempo_mn);
        [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm3, epsilon, count_success_newton, count_failure_newton);

        % 2) Modified Newton + Prec
        tic;
        [xk3_prec, fk3_prec, gradfk_norm3_prec, k3_prec, xseq3_prec, btseq3_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad , problem_213_hess, ...
            kmax, tolgrad, c1, rho, btmax);
        tempo_mn_prec = toc;
        rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton + Prec", k3_prec, norm(xk3_prec), fk3_prec, gradfk_norm3_prec, tempo_mn_prec);
        [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm3_prec, epsilon, count_success_newton, count_failure_newton);

        
        % 3) Modified Newton FD
        tic;
        [xk3_fd, fk3_fd, gradfk_norm3_fd, k3_fd, xseq3_fd, btseq3_fd] = ...
            modified_newton_backtracking(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad_fd , problem_213_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd = toc;
        rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton FD", k3_fd, norm(xk3_fd), fk3_fd, gradfk_norm3_fd, tempo_mn_fd);
        [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm3_fd, epsilon, count_success_newton, count_failure_newton);
    
        % 4) Modified Newton FD + Prec
        tic;
        [xk3_fd_prec, fk3_fd_prec, gradfk_norm3_fd_prec, k3_fd_prec, xseq3_fd_prec, btseq3_fd_prec] = ...
            modified_newton_backtracking_preconditioning(x_bar_problem_213, problem_213_fun, ...
            problem_213_grad_fd , problem_213_hess_fd, ...
            kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
        tempo_mn_fd_prec = toc;
        rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton FD + Prec", k3_fd_prec, norm(xk3_fd_prec), fk3_fd_prec, gradfk_norm3_fd_prec, tempo_mn_fd_prec);
        [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm3_fd_prec, epsilon, count_success_newton, count_failure_newton);
        
        
        % With 10 starting points generated with uniform distribution in a hyper-cube
        for i = 1:num_points
            
            x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;
            startLabel = sprintf('%d', i); % tabella: 1..10

            % Modified Newton
            tic;
            [xk_rand, fk_rand, gradfk_norm_rand, k_rand] = ...
                modified_newton_backtracking(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);
            tempo_mn_rand = toc;
            rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton", k_rand, norm(xk_rand), fk_rand, gradfk_norm_rand, tempo_mn_rand);
            [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm_rand, epsilon, count_success_newton, count_failure_newton);

            % Newton with preconditioning
            tic;
            [xk_rand_prec, fk_rand_prec, gradfk_norm_rand_prec, k_rand_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad , problem_213_hess, ...
                kmax, tolgrad, c1, rho, btmax);
            tempo_mn_prec_rand = toc;
            rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton + Prec", k_rand_prec, norm(xk_rand_prec), fk_rand_prec, gradfk_norm_rand_prec, tempo_mn_prec_rand);
            [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm_rand_prec, epsilon, count_success_newton, count_failure_newton);
            
            % Finite difference
            tic;
            [xk_rand_fd, fk_rand_fd, gradfk_norm_rand_fd, k_rand_fd, xseq_fd, btseq_fd] = ...
                modified_newton_backtracking(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            tempo_mn_fd_rand = toc;
            rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton FD", k_rand_fd, norm(xk_rand_fd), fk_rand_fd, gradfk_norm_rand_fd, tempo_mn_fd_rand);
            [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm_rand_fd, epsilon, count_success_newton, count_failure_newton);

            % Finite difference + prec
            tic;
            [xk_rand_fd_prec, fk_rand_fd_prec, gradfk_norm_rand_fd_prec, k_rand_fd_prec, xseq_fd_prec, btseq_fd_prec] = ...
                modified_newton_backtracking_preconditioning(x0_i, problem_213_fun, ...
                problem_213_grad_fd , problem_213_hess_fd, ...
                kmax, tolgrad, c1, rho, btmax, hstep, hstep_i);
            tempo_mn_fd_prec_rand = toc;
            rows(end+1) = experiment_utils('make_row',startLabel, "Modified Newton FD + Prec", k_rand_fd_prec, norm(xk_rand_fd_prec), fk_rand_fd_prec, gradfk_norm_rand_fd_prec, tempo_mn_fd_prec_rand);
            [count_success_newton, count_failure_newton] = experiment_utils('update_counts',gradfk_norm_rand_fd_prec, epsilon, count_success_newton, count_failure_newton);

        end

        % ------------------ Print the table -------------------------
        experiment_utils('print_table', fid, n, hstep, rows);
    end

     % Summary counts (per n)
    fprintf(fid, "\nSUMMARY for n = %d\n", n);
    fprintf(fid, "Successes (||grad|| < %.1e): %d\n", epsilon, count_success_newton);
    fprintf(fid, "Failures  (||grad|| >= %.1e): %d\n\n", epsilon, count_failure_newton);

end


% ======================= NELDER-MEAD ===========================

fprintf(fid, "\n\n============================================================\n");
fprintf(fid, "NELDER-MEAD method - Tables\n");
fprintf(fid, "============================================================\n\n");

for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;

    fprintf('Running Nelder-Mead for n = %d\n', n);

    % First starting point x_bar
    x_bar_problem_213 = ones(n,1);

    % Prepare table rows for THIS n
    % Columns: start | method | iter | xnorm | fval | D | time
    rows_nm = struct('start',{},'method',{},'iter',{},'xnorm',{},'fval',{},'D',{},'time',{});

    % ------------------ Starting point: x_bar ------------------
    startLabel = 'x_bar';

    tic;
    simplex_problem_213 = nelder_mead_n(x_bar_problem_213, problem_213_fun, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
    tempo_nelder_mead = toc; 
    D = max(pdist(simplex_problem_213'));   
    x_best = simplex_problem_213(:,1);
    rows_nm(end+1) = experiment_utils('make_row_nm', startLabel, "Nelder-Mead", NaN, norm(x_best), ...
                                  problem_213_fun(x_best), D, tempo_nelder_mead);
    [count_success_nelder, count_failure_nelder] = ...
    experiment_utils('update_counts', D, epsilon, count_success_nelder, count_failure_nelder);


    

    % With 10 starting points generated with uniform distribution in a hyper-cube
    for i = 1:num_points
        x0_i = x_bar_problem_213 + 2 * rand(n,1) - 1;

        tic;
        simplex_i = nelder_mead_n(x0_i, problem_213_fun, n , ...
            rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
        tempo_nelder_mead = toc; 
        Di = max(pdist(simplex_i'));
        x_best_i = simplex_i(:,1);
        rows_nm(end+1) = experiment_utils('make_row_nm', startLabel, "Nelder-Mead", NaN, norm(x_best_i), ...
                                  problem_213_fun(x_best_i), D, tempo_nelder_mead);
        [count_success_nelder, count_failure_nelder] = ...
        experiment_utils('update_counts', D, epsilon, count_success_nelder, count_failure_nelder);
       
    end
    experiment_utils('print_table_nelder', fid, n, rows_nm);

    % Summary counts (per n)
    fprintf(fid, "\nSUMMARY for n = %d\n", n);
    fprintf(fid, "Successes (||grad|| < %.1e): %d\n", epsilon, count_success_nelder);
    fprintf(fid, "Failures  (||grad|| >= %.1e): %d\n\n", epsilon, count_failure_nelder);
    
end

fclose(fid);



