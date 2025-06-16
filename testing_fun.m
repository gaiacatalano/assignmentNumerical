clc
clear
close all

%Funzione di Rosenbrock
f_rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;

% Definizione della griglia
[x, y] = meshgrid(-2:0.05:2, -1:0.05:3);
z = 100*(y - x.^2).^2 + (1 - x).^2;

% Grafico 3D
figure;
surf(x, y, z);
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('Funzione di Rosenbrock - Surface 3D');
shading interp;
colormap turbo;
view(45, 30);
grid on;

% Starting points
x0_a = [1.2; 1.2];
x0_b = [-1.2; 1];

% Backtracking parameters
rho = 0.5;

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
tempo_nm_a = toc;

tic;
[simplex_nm_b, k_nm_b] = nelder_mead(x0_b, f_rosenbrock, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol);
tempo_nm_b = toc;


% MODIFIED NEWTON

tolgrad = 1e-8;
gradf = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                200*(x(2) - x(1)^2) ];
Hessf = @(x) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
               -400*x(1), 200 ];
c1 = 1e-4;
btmax = 20;

tic;
[xk_a, fk_a, gradfk_norm_a, k_a, xseq_a, btseq_a] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
tempo_mn_a = toc;
x_newton_a = xk_a;

tic;
[xk_b, fk_b, gradfk_norm_b, k_b, xseq_b, btseq_b] = ...
    modified_newton_bcktrck(x0_a, f_rosenbrock, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax);
tempo_mn_b = toc;
x_newton_b = xk_b;


fprintf(fid, "== RISULTATI PER %s ==\n", "MODIFIED NEWTON WITH BACKTRACKING");
fprintf(fid, '\n');

fprintf(fid, "Punto di partenza: [%f, %f]\n", x0_a(1), x0_a(2));
fprintf(fid, "Numero di iterazioni: %d\n", k_a);
fprintf(fid, "Valore funzione finale: %.10f\n", f_rosenbrock(xk_a));
fprintf(fid, "Tempo di esecuzione: %.6f secondi\n", tempo_mn_a);
fprintf(fid, "Norma gradiente finale: %.2e\n", gradfk_norm_a);
fprintf(fid, "Numero totale di backtracking: %d\n", sum(btseq_a));
fprintf(fid, '\n');


fprintf(fid, "Punto di partenza: [%f, %f]\n", x0_b(1), x0_b(2));
fprintf(fid, "Numero di iterazioni: %d\n", k_b);
fprintf(fid, "Valore funzione finale: %.10f\n", f_rosenbrock(xk_b));
fprintf(fid, "Tempo di esecuzione: %.6f secondi\n", tempo_mn_b);
fprintf(fid, "Norma gradiente finale: %.2e\n", gradfk_norm_b);
fprintf(fid, "Numero totale di backtracking: %d\n", sum(btseq_b));

fprintf(fid, '\n');
fprintf(fid, "\n== RISULTATI PER NELDER-MEAD ==\n");
fprintf(fid, '\n');

fprintf(fid, "Punto di partenza: [%f, %f]\n", x0_a(1), x0_a(2));
fprintf(fid, "Numero di iterazioni: %d\n", k_nm_a);
fprintf(fid, "Valore funzione finale: %.10f\n", f_rosenbrock(simplex_nm_a(:,1)));
fprintf(fid, "Tempo di esecuzione: %.6f secondi\n", tempo_nm_a);

fprintf(fid, '\n');
fprintf(fid, "Punto di partenza: [%f, %f]\n", x0_b(1), x0_b(2));
fprintf(fid, "Numero di iterazioni: %d\n", k_nm_b);
fprintf(fid, "Valore funzione finale: %.10f\n", f_rosenbrock(simplex_nm_b(:,1)));
fprintf(fid, "Tempo di esecuzione: %.6f secondi\n", tempo_nm_b);

methods = {'Newton x0_a', 'Newton x0_b', 'Nelder-Mead x0_a', 'Nelder-Mead x0_b'};

times = [tempo_mn_a, tempo_mn_b, tempo_nm_a, tempo_nm_b];

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