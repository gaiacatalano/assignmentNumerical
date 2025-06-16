% Intervallo di valori per x1 e x2
x1 = linspace(-2, 2, 100);
x2 = linspace(-1, 3, 100);

% Creazione della griglia
[X1, X2] = meshgrid(x1, x2);

% Definizione della funzione di Rosenbrock
F = 100 * (X1.^2 - X2).^2 + (X1 - 1).^2;

% Plot della superficie
figure;
surf(X1, X2, F)
xlabel('x_1')
ylabel('x_2')
zlabel('F(x_1, x_2)')
title('Chained Rosenbrock Function (2D)')
shading interp
colormap turbo
colorbar
