clc
clear
close all

%Funzione di Rosenbrock
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

% Nelder-Mead parameters
rho_nm = 1;
chi_nm = 2;
gamma_nm = 0.5;
sigma_nm = 0.5;

fid = fopen('output_rosenbrock.txt', 'w');

tic;
simplex_nm = nelder_mead(x0_a, f_rosenbrock, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
tempo_nm = toc;
fprintf(fid, "Tempo di esecuzione Nelder mead %d\n", tempo_nm);

% Newton parameters
tolgrad = 1e-8;
gradf = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                200*(x(2) - x(1)^2) ];
Hessf = @(x) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
               -400*x(1), 200 ];
c1 = 1e-8;
btmax = 20;

tic;
[xk, fk, gradfk_norm, k, xseq, btseq] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
tempo_mn = toc;
fprintf(fid, "Tempo di esecuzione Modified Newton %d\n", tempo_mn);
x_newton = xk;

% Anche cambiando tolleranza e numero di iterazioni, il simplesso rimane lo
% stesso; idem con il metodo di Newton

fclose(fid);