function g = chained_rosenbrock_grad_fd(x, h, bool_hstep_i)

    n = length(x);
    g = zeros(n, 1);

    if bool_hstep_i == 1
        hstep_i = abs(x)*h;
    else
        hstep_i = ones(n,1) * h;
    end

    % Funzione locale f_i(x): solo i=2 to n
    for i = 1:n
        %fx = 0;
        %xi=x(i);

        hi = hstep_i(i);
        e = zeros(n,1);
        e(i) = 1;

        fp = chained_rosenbrock_fvalue(x + hi*e);
        fm = chained_rosenbrock_fvalue(x - hi*e);

        g(i) = (fp - fm) / (2*hi);
        % evito di calcolare tutt F(x) n volte ma calcolo solo i due valori
        % f_i che mi servono

        %calcolo di f_i
        % if i >= 2
        %     xim1 = x(i-1);
        %     g(i) = g(i) + 100*(h-2*(xim1^2 - xi));  %già tutto diviso per h
        % 
        % end
        % 
        % %calcolo di f_{i+1}
        % if i <= n - 1
        %     xip1   = x(i+1);
        %     g(i) = g(i) + 100*(h^3 + 4*xi^2*h + 4*xi*h^2 + 2*(h+2*xi)*(xi^2-xip1)) +h +2*(xi-1); %già tutto diviso per h
        % end
        
        % if i == 1
        %     g(i) = 400*x(i)*(h^2 + x(i)^2 - x(i+1)) + 2*(x(i) - 1);
        % elseif i == n
        %     g(i) = -200*(x(i-1)^2 - x(i));
        % else
        %     g(i) = 400*x(i)*(h^2 + x(i)^2 - x(i+1)) - 200*(x(i-1)^2 - x(i)) + 2*(x(i) - 1);
        % end
        

    end
end