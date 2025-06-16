% Intervallo per x1 e x2
[x1, x2] = meshgrid(linspace(-2, 2, 200), linspace(-1, 3, 200));

% Funzione di Rosenbrock
f = 100 * (x2 - x1.^2).^2 + (1 - x1).^2;

% Plot 3D
figure
surf(x1, x2, f)
xlabel('x_1')
ylabel('x_2')
zlabel('f(x_1, x_2)')
title('Funzione di Rosenbrock')
shading interp
