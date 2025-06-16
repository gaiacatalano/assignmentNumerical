clc
clear
close all

f_rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Plot of the function and the minimun point
[x, y] = meshgrid(-2:0.05:2, -1:0.05:3);
z = 100*(y - x.^2).^2 + (1 - x).^2;

figure;
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Rosenbrock function and its minimun point');
shading interp;
colormap turbo;
view(45, 30);
hold on;
plot3(1, 1, f_rosenbrock([1; 1]), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Starting points used for the test
x0_a = [1.2; 1.2];
x0_b = [-1.2; 1];

% Stopping parameters
tol = 1e-8;
kmax = 2000;

% Nelder-Mead parameters
rho_nm = 1;
chi_nm = 2;
gamma_nm = 0.5;
sigma_nm = 0.5;

fid = fopen('output_rosenbrock.txt', 'w');

% NELDER-MEAD

tic;
[simplex_nm_a, k_nm_a] = nelder_mead(x0_a, f_rosenbrock, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
time_nm_a = toc;

tic;
[simplex_nm_b, k_nm_b] = nelder_mead(x0_b, f_rosenbrock, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
time_nm_b = toc;


% MODIFIED NEWTON

% Backtracking parameters
rho = 0.5;
tolgrad = 1e-8;
c1 = 1e-4;
btmax = 20;

gradf = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                200*(x(2) - x(1)^2) ];
Hessf = @(x) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
               -400*x(1), 200 ];

tic;
[xk_a, fk_a, gradfk_norm_a, k_a, xseq_a, btseq_a] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
time_mn_a = toc;
x_newton_a = xk_a;

tic;
[xk_b, fk_b, gradfk_norm_b, k_b, xseq_b, btseq_b] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
time_mn_b = toc;
x_newton_b = xk_b;


fprintf(fid, "== RESULTS FOR %s ==\n", "MODIFIED NEWTON WITH BACKTRACKING");
fprintf(fid, '\n');

fprintf(fid, "Starting point: [%f, %f]\n", x0_a(1), x0_a(2));
fprintf(fid, "Number of iterations: %d\n", k_a);
fprintf(fid, "Final value of the function: %.10f\n", f_rosenbrock(xk_a));
fprintf(fid, "Execution time: %.6f secondi\n", time_mn_a);
fprintf(fid, "Final norm of the gradient: %.2e\n", gradfk_norm_a);
fprintf(fid, "Total number of backtracking: %d\n", sum(btseq_a));
fprintf(fid, '\n');


fprintf(fid, "Starting point: [%f, %f]\n", x0_b(1), x0_b(2));
fprintf(fid, "Number of iterations: %d\n", k_b);
fprintf(fid, "Final value of the function: %.10f\n", f_rosenbrock(xk_b));
fprintf(fid, "Execution time: %.6f secondi\n", time_mn_b);
fprintf(fid, "Final norm of the gradient: %.2e\n", gradfk_norm_b);
fprintf(fid, "Total number of backtracking: %d\n", sum(btseq_b));

fprintf(fid, '\n');
fprintf(fid, "\n== RISULTS FOR NELDER-MEAD ==\n");
fprintf(fid, '\n');

fprintf(fid, "Starting point: [%f, %f]\n", x0_a(1), x0_a(2));
fprintf(fid, "Number of iterations: %d\n", k_nm_a);
fprintf(fid, "Final value of the function: %.10f\n", f_rosenbrock(simplex_nm_a(:,1)));
fprintf(fid, "Execution time: %.6f secondi\n", time_nm_a);

fprintf(fid, '\n');
fprintf(fid, "Starting point: [%f, %f]\n", x0_b(1), x0_b(2));
fprintf(fid, "Number of iterations: %d\n", k_nm_b);
fprintf(fid, "Final value of the function: %.10f\n", f_rosenbrock(simplex_nm_b(:,1)));
fprintf(fid, "Execution time: %.6f secondi\n", time_nm_b);

methods = {'Newton x0_a', 'Newton x0_b', 'Nelder-Mead x0_a', 'Nelder-Mead x0_b'};

times = [time_mn_a, time_mn_b, time_nm_a, time_nm_b];

figure;
bar(times);
set(gca, 'XTickLabel', methods);
ylabel('Tempo di esecuzione (s)');
title('Confronto dei tempi di esecuzione');
grid on;

iterations = [k_a, k_b, k_nm_a, k_nm_b];  % Nelder-Mead ha fisso kmax nel tuo caso

figure;
bar(iterations);
set(gca, 'XTickLabel', methods);
ylabel('Numero di iterazioni');
title('Confronto del numero di iterazioni');
grid on;




fclose(fid);