function g = discrete_boundary_value_grad_fd(x, hstep)
    if nargin < 2
        hstep = 1e-6;
    end

    n = length(x);
    h = 1 / (n + 1);  
    g = zeros(n,1);
    x_ext = [0; x; 0];  % estensione con x0 = 0, x_{n+1} = 0

    for i = 1:n
        xim1 = x_ext(i);
        xi = x_ext(i+1);
        xip1 = x_ext(i+2);
        fx=0;

        if i > 2
            xim2 = x_ext(i-1);
            fim1 = 2*xim1 - xim2 - xi + (h^2 / 2)*(xim1 + (i-1)*h + 1)^3;
            %fx = fx + fi^2;
            g(i) = g(i) - 2*fim1 + hstep;
        end
        
        fi = 2*xi - xim1 - xip1 + (h^2 / 2)*(xi + (i)*h + 1)^3;
        %fx = fx + fi^2;
        g(i) = g(i) + (2+((h*hstep)^2)/2 + ((3*h^2)/2)*((xi+i*h+1)^2 + (xi+i*h+1)*hstep))^2 + 2*fi*(2+((h*hstep)^2)/2 + ((3*h^2)/2)*((xi+i*h+1)^2 + (xi+i*h+1)*hstep));

        if i < n
            xip2 = x_ext(i+3);
            fip1 = 2*xip1 - xi - xip2 + (h^2 / 2)*(xip1 + (i+1)*h + 1)^3;
            %fx = fx + fi^2;
            g(i) = g(i) + hstep - 2*fip1;
        end

    end
end