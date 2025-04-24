clc
clear
close all

% Seed
rng(349131);

% Dimension
d = 3:1:5;
% n = [10,25,50];
n = 10;

num_points = 10;

% Problem 32
x_bar = -1;
starting_points = 2*rand(n, num_points) - 1;  % genera valori in [-1, 1]
starting_points = starting_points + x_bar;    % trasla il centro in x

% Funzione
F = broyden_tridiagonal(n, x_bar);

% Nelder-Mead

% Stopping parameters
tol = 1e-8;
kmax = 2000;

% Nelder-Mead parameters
rho_nm = 1;
chi_nm = 2;
gamma_nm = 0.5;
sigma_nm = 0.5;

simplex_nm = nelder_mead_n(x_bar*ones(n,1), F, n , rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);

% Modified

% Backtracking parameters
rho = 0.5;
c = 10e-4;

% Newton parameters
tolgrad = 1e-8;
gradf = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                200*(x(2) - x(1)^2) ];
Hessf = @(x) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
               -400*x(1), 200 ];
c1 = 1e-8;
btmax = 20;

[xk, fk, gradfk_norm, k, xseq, btseq] = ...
    modified_newton_bcktrck(x_bar*ones(1,n), F, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
x_newton = xk
