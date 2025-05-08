function g = chained_rosenbrock_grad_fd(x, h)
    if nargin < 2
        h = 1e-6;
    end

    n = length(x);
    g = zeros(n, 1);

    % Funzione locale f_i(x): solo i=2 to n
    for i = 1:n
        fx = 0;

        % evito di calcolare tutt F(x) n volte ma calcolo solo i due valori
        % f_i che mi servono

        %calcolo di f_i
        if i >= 2
            xim1 = x(i-1);
            xi   = x(i);
            f_i = 100 * (xim1^2 - xi)^2 + (xim1 - 1)^2;
            fx = fx + f_i;
        end

        %calcolo di f_{j+1}
        if i <= n - 1
            xim1 = x(i);
            xi   = x(i+1);
            f_i = 100 * (xim1^2 - xi)^2 + (xim1 - 1)^2;
            fx = fx + f_i;
        end

        % Punto perturbato
        xh = x;
        xh(i) = xh(i) + h;

        fxh = 0;
        if i >= 2
            xim1 = xh(i-1);
            xi   = xh(i);
            f_i = 100 * (xim1^2 - xi)^2 + (xim1 - 1)^2;
            fxh = fxh + f_i;
        end
        if i <= n - 1
            xim1 = xh(i);
            xi   = xh(i+1);
            f_i = 100 * (xim1^2 - xi)^2 + (xim1 - 1)^2;
            fxh = fxh + f_i;
        end

        g(i) = (fxh - fx) / h;
    end
end