% Parameters
n = 2;
h = 1 / (n + 1);

% Grid for x1 e x2
x1 = linspace(-1, 1, 100);
x2 = linspace(-1, 1, 100);
[X1, X2] = meshgrid(x1, x2);

x0 = 0;
x3 = 0;

F = (2*X1 - x0 - X2 + h^2 * (X1 + 1*h + 1).^(3/2)).^2 + ...
    (2*X2 - X1 - x3 + h^2 * (X2 + 2*h + 1).^(3/2)).^2;

% Plot
figure;
surf(X1, X2, F)
xlabel('x_1')
ylabel('x_2')
zlabel('F(x)')
title('Discrete Boundary Value Problem (n = 2)')
shading interp
colormap parula
colorbar