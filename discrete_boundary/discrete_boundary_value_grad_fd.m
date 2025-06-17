% Approximates the gradient of the Discrete Boundary Value problem 
% using finite difference methods

function g = discrete_boundary_value_grad_fd(x, hstep, bool_hstep_i)


    if nargin < 2
        hstep = 1e-6;
    end

    if bool_hstep_i==1
        hstep_i= abs(x)*hstep;
    end

    n = length(x);
    h = 1 / (n + 1);  
    g = zeros(n,1);
    x_ext = [0; x; 0]; 

    for i = 1:n
        xim1 = x_ext(i);
        xi = x_ext(i+1);
        xip1 = x_ext(i+2);

        if bool_hstep_i==1
            hstep = hstep_i(i);
        end

        if i > 1
            xim2 = x_ext(i-1);
            fim1 = 2*xim1 - xim2 - xi + (h^2 / 2)*(xim1 + (i-1)*h + 1)^3;

            if bool_hstep_i==0
                g(i) = g(i) - 2*fim1 + hstep;
            else
                g(i) = g(i) + hstep_i(i-1) - 2*fim1;
            end
        end
        
        fi = 2*xi - xim1 - xip1 + (h^2 / 2)*(xi + (i)*h + 1)^3;
        g(i) = g(i) + (2+((h*hstep)^2)/2 + ((3*h^2)/2)*((xi+i*h+1)^2 + ...
            (xi+i*h+1)*hstep))^2 + 2*fi*(2+((h*hstep)^2)/2 + ...
            ((3*h^2)/2)*((xi+i*h+1)^2 + (xi+i*h+1)*hstep));

        if i < n
            xip2 = x_ext(i+3);
            fip1 = 2*xip1 - xi - xip2 + (h^2 / 2)*(xip1 + (i+1)*h + 1)^3;

            if bool_hstep_i==0
                g(i) = g(i) + hstep - 2*fip1;
            else
                g(i) = g(i) + hstep_i(i+1) - 2*fip1;
            end
        end

    end
end