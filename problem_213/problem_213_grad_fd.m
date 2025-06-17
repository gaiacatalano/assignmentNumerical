% Approximates the gradient of the Problem 213
% using finite difference methods

function g = problem_213_grad_fd(x, hstep, bool_hstep_i)

    if nargin < 2
        hstep = 1e-6;
    end

    if bool_hstep_i == 1
        hstep_i = abs(x)*hstep;
    end

    n = length(x);
    h = 1 / (n + 1);  
    g = zeros(n,1);
    x_ext = [0; x; 0]; 

    for i = 1:n
        xim1 = x_ext(i);
        xi = x_ext(i+1);
        xip1 = x_ext(i+2);
        fx=0;

        if bool_hstep_i==1
            hstep = hstep_i(i);
        end

        if i > 2
            xim2 = x_ext(i-1);
            fim1 = 2*xim1 + (h^2)*(xim1 + sin(xim1))-xim2 - xi;
            g(i) = g(i) + 0.5*hstep - fim1;
        end

        fi = 2*xi + (h^2)*(xi + sin(xi))-xim1 - xip1;
        g(i) = g(i) + 0.5*(4*hstep + (h^4)*hstep + 4*(h^2)*sin(xi) + 4*hstep*h^2 + 4*(h^2)*sin(xi+hstep) - 2*(h^4)*sin(xi) + 2*(h^4)*sin(xi+hstep)) + (1/(2*hstep))*(h^4)*((sin(xi)^2) + (sin(xi+hstep)^2) - 2*sin(xi)*sin(xi+hstep)) + fi*(2+h^2) + (fi*(h^2)/hstep)*(sin(xi+hstep) - sin(xi));

        if i < n
            xip2 = x_ext(i+3);
            fip1 = 2*xip1 + (h^2)*(xip1 + sin(xip1))-xi - xip2;
            g(i) = g(i) + 0.5*hstep - fip1;
        end
    end
end
