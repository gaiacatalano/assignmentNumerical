%Funzione di Rosenbrock
%f_rosenbrock = @(x1, x2) 100*(x2-x1^2)^2 + (1-x1)^2;
f_rosenbrock = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;


%Starting points
x0_a = [1.2; 1.2];
x0_b = [-1.2, 1];

%Backtracking parameters
rho = 0.5;
c = 10e-4;

%Stopping parameters
tol_x = 10e-5;
kmax = 100;
max_no_improvement = 10;

%Nelder-Mead parameters
rho_nm = 1;
chi_nm = 2;
gamma_nm = 0.5;
sigma_nm = 0.5;

NM = nelder_mead(x0_a, f_rosenbrock, 2, rho_nm, chi_nm, gamma_nm, sigma_nm, kmax, tol_x, max_no_improvement)
