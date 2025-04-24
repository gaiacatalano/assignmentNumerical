clc
clear
close all

%Funzione di Rosenbrock
%f_rosenbrock = @(x1, x2) 100*(x2-x1^2)^2 + (1-x1)^2;
f_rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Starting points
x0_a = [1.2; 1.2];
x0_b = [-1.2; 1];

% Backtracking parameters
rho = 0.5;
c = 10e-4;

% Stopping parameters
tol = 1e-8;
kmax = 2000;
max_no_improvement = 40;

% Nelder-Mead parameters
rho_nm = 1;
chi_nm = 2;
gamma_nm = 0.5;
sigma_nm = 0.5;
tol_f = 1e-16;

simplex_nm = nelder_mead(x0_a, f_rosenbrock, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol)

% Newton parameters
tolgrad = 1e-8;
gradf = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                200*(x(2) - x(1)^2) ];
Hessf = @(x) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
               -400*x(1), 200 ];
c1 = 1e-8;
btmax = 20;

[xk, fk, gradfk_norm, k, xseq, btseq] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
x_newton = xk

% Anche cambiando tolleranza e numero di iterazioni, il simplesso rimane lo
% stesso; idem con il metodo di Newton
