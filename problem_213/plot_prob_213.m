n = 10;
h = 1 / (n + 1);

% Grid for x4 e x5
[x4_vals, x5_vals] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
F_vals = zeros(size(x4_vals));

for i = 1:size(x4_vals, 1)
    for j = 1:size(x4_vals, 2)
        x = ones(n, 1);
        x(4) = x4_vals(i, j);
        x(5) = x5_vals(i, j);
        
        x_ext = [0; x; 1];
        
        F = 0;
        for k = 1:n
            f_k = 2*x_ext(k+1) + h^2 * (x_ext(k+1) + sin(x_ext(k+1))) ...
                  - x_ext(k) - x_ext(k+2);
            F = F + f_k^2;
        end
        F_vals(i, j) = 0.5 * F;
    end
end

% Plot 
surf(x4_vals, x5_vals, F_vals)
xlabel('x_4')
ylabel('x_5')
zlabel('F(x)')
title('Plot di F(x) in funzione di x_4 e x_5')
shading interp
