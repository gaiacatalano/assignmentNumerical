function g = chained_rosenbrock_grad_fd(x, h, bool_hstep_i)
    if nargin < 2
        h = 1e-6;
    end

    if bool_hstep_i == 1
        hstep_i = abs(x)*h;
    end

    n = length(x);
    g = zeros(n, 1);

    % Funzione locale f_i(x): solo i=2 to n
    for i = 1:n
        %fx = 0;
        xi=x(i);

        if bool_hstep_i == 1
            h = hstep_i(i);
        end

        % evito di calcolare tutt F(x) n volte ma calcolo solo i due valori
        % f_i che mi servono

        %calcolo di f_i
        if i >= 2
            xim1 = x(i-1);
            g(i) = g(i) + 100*(h-2*(xim1^2 - xi));  %già tutto diviso per h
            
        end

        %calcolo di f_{i+1}
        if i <= n - 1
            xip1   = x(i+1);
            g(i) = g(i) + 100*(h^3 + 4*xi^2*h + 4*xi*h^2 + 2*(h+2*xi)*(xi^2-xip1)) +h +2*(xi-1); %già tutto diviso per h
        end

        

    end
end