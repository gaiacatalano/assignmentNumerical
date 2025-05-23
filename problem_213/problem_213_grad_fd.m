function g = problem_213_grad_fd(x, hstep)
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
    %     if i > 1
    %         xim1 = x_ext(i-1);
    %         xi   = x_ext(i);
    %         xip1 = x_ext(i+1);
    %         fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %         fx = fx + fi^2;            
    %     end
    % 
    %     xim1 = x_ext(i);
    %     xi   = x_ext(i+1);
    %     xip1 = x_ext(i+2);
    %     fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %     fx = fx + fi^2;
    % 
    %     if i < n
    %         xim1 = x_ext(i+1);
    %         xi   = x_ext(i+2);
    %         xip1 = x_ext(i+3);
    %         fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %         fx = fx + fi^2;
    %     end
    % 
    %     %pezzo perturbato
    %     xh = x;
    %     xh(i) = xh(i) + hstep;
    %     xh_ext = [0; xh; 0];
    % 
    %     fxh = 0;
    %     if i > 1
    %         xim1 = xh_ext(i-1);
    %         xi   = xh_ext(i);
    %         xip1 = xh_ext(i+1);
    %         fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %         fxh = fxh + fi^2;            
    %     end
    % 
    %     xim1 = xh_ext(i);
    %     xi   = xh_ext(i+1);
    %     xip1 = xh_ext(i+2);
    %     fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %     fxh = fxh + fi^2;
    % 
    %     if i < n
    %         xim1 = xh_ext(i+1);
    %         xi   = xh_ext(i+2);
    %         xip1 = xh_ext(i+3);
    %         fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
    %         fxh = fxh + fi^2;
    %     end
    % 
    %     g(i) = (fxh - fx) / hstep;
    % 
    % end


end
