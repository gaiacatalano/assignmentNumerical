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
kmax = 2000;
    
% Ciclo for per valori di n adatti a Newton
for p=1:length(d)

    n = 10^d(p);

    % Backtracking parameters
    rho = 0.5;
    c = 10e-4;
    
    % Newton parameters
    tolgrad = 1e-5;
    c1 = 1e-8;
    btmax = 20;


    % Problem 31
    x_bar = -1;

    % % Starting points generati random
    % starting_points = 2*rand(n, num_points) - 1;  % genera valori in [-1, 1]
    % starting_points = starting_points + x_bar;    % trasla il centro in x
    
    [F, grad, H] = broyden_tridiagonal(n, x_bar);    
    [xk, fk, gradfk_norm, k, xseq, btseq] = ...
        modified_newton_bcktrck(x_bar*ones(n,1), F, grad , H, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton = xk;


    % Problem 14
    x_bar2 = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar2(i) = i*h*(1-i*h);
    end
    [F2, grad2, H2] = discrete_boundary_value_problem(n,x_bar2);
    [xk2, fk2, gradfk_norm2, k2, xseq2, btseq2] = ...
        modified_newton_bcktrck(x_bar2, F2, grad2 , H2, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton2 = xk2;

    % Problem 1
    x_bar3 = zeros(n,1);         
    x_bar3(mod(1:n,2)==1) = -1.2; 
    x_bar3(mod(1:n,2)==0) = 1.0;  
    [F3, grad3, H3] = chained_rosenbrock(n,x_bar3);
    [xk3, fk3, gradfk_norm3, k3, xseq3, btseq3] = ...
        modified_newton_bcktrck(x_bar3, F3, grad3 , H3, ...
        kmax, tolgrad, c1, rho, btmax);
    x_newton3 = xk3;

    % Restituisco valori newton dei problemi
    x_newton
    x_newton2
    x_newton3
end

% Ciclo for per valori di n adatti a Nelder-Mead
for n = [10,25,50]

    % Nelder-Mead parameters
    rho_nm = 1;
    chi_nm = 2;
    gamma_nm = 0.5;
    sigma_nm = 0.5;
    
 
    % Problem 31
    x_bar = -1;

    % % Starting points generati random
    % starting_points = 2*rand(n, num_points) - 1;  % genera valori in [-1, 1]
    % starting_points = starting_points + x_bar;    % trasla il centro in x

    [F, grad, H] = broyden_tridiagonal(n, x_bar);
    simplex_nm = nelder_mead_n(x_bar*ones(n,1), F, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
       

    % Problem 14
    x_bar2 = zeros(n,1);
    h = 1/(n+1);
    for i=1:n
        x_bar2(i) = i*h*(1-i*h);
    end
    [F2, grad2, H2] = discrete_boundary_value_problem(n,x_bar2);
    simplex_nm2 = nelder_mead_n(x_bar2, F2, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

    % Problema 1
    x_bar3 = zeros(n,1);        
    x_bar3(mod(1:n,2)==1) = -1.2; 
    x_bar3(mod(1:n,2)==0) = 1.0; 
    [F3, grad3, H3] = chained_rosenbrock(n,x_bar3);
    simplex_nm3 = nelder_mead_n(x_bar3, F3, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);


    % Restituisco valore migliore di ogni simplesso
    simplex_nm(:,1)
    simplex_nm2(:,1)
    simplex_nm3(:,1)
end