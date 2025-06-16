% Parametri
n = 2;
h = 1 / (n + 1);

% Definizione della griglia per x1 e x2
x1 = linspace(-1, 1, 100);
x2 = linspace(-1, 1, 100);
[X1, X2] = meshgrid(x1, x2);

% Condizioni al bordo
x0 = 0;
x3 = 0;

% Calcolo del valore di F(x) per n=2:
% i = 1 --> usa x0, x1, x2
% i = 2 --> usa x1, x2, x3

F = (2*X1 - x0 - X2 + h^2 * (X1 + 1*h + 1).^(3/2)).^2 + ...
    (2*X2 - X1 - x3 + h^2 * (X2 + 2*h + 1).^(3/2)).^2;

% Plot della superficie
figure;
surf(X1, X2, F)
xlabel('x_1')
ylabel('x_2')
zlabel('F(x)')
title('Discrete Boundary Value Problem (n = 2)')
shading interp
colormap parula
colorbar
