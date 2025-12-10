% Approximates the gradient of the Problem 213
% using finite difference methods

function g = problem_213_grad_fd(x, hstep, bool_hstep_i)

    n = length(x);
    if nargin < 2 || isempty(hstep)
        hstep = 1e-6;
    end

    if bool_hstep_i == 1
        hstep_i = abs(x)*hstep;
   else
        % step fisso
        hstep_i = hstep * ones(n,1);
   
    end

    g = zeros(n,1);

    %h = 1 / (n + 1);  
    %x = [0; x; 1]; 

    for i = 1:n
        hi = hstep_i(i);

        e = zeros(n,1);
        e(i) = 1;

        Fp  = problem_213_fvalue(x + hi*e);
        Fm = problem_213_fvalue(x - hi*e);

        g(i) = (Fp - Fm) / (2*hi);
    end
end


%         xim1 = x(i);
%         xi = x(i+1);
%         xip1 = x(i+2);
%         %fx=0;
% 
%         if bool_hstep_i==1
%             hstep = hstep_i(i);
%         end
% 
%         if i > 1
%             xim2 = x(i-1);
%             fim1 = 2*xim1 + (h^2)*(xim1 + sin(xim1))-xim2 - xi;
%             %g(i) = g(i) + 0.5*hstep - fim1;
%             g(i) = g(i) - fim1;
%         end
% 
%         fi = 2*xi + (h^2)*(xi + sin(xi))-xim1 - xip1;
%         %g(i) = g(i) + 0.5*(4*hstep + (h^4)*hstep + 4*(h^2)*sin(xi) + 4*hstep*h^2 + 4*(h^2)*sin(xi+hstep) - 2*(h^4)*sin(xi) + 2*(h^4)*sin(xi+hstep)) + (1/(2*hstep))*(h^4)*((sin(xi)^2) + (sin(xi+hstep)^2) - 2*sin(xi)*sin(xi+hstep)) + fi*(2+h^2) + (fi*(h^2)/hstep)*(sin(xi+hstep) - sin(xi));
%         g(i) = g(i) + 1/(4*hstep)*(h^4*(sin(x(i+1)+hstep)^2 - sin(x(i+1)-hstep)^2 - 2*sin(x(i+1))*(sin(x(i+1)+hstep) - sin(x(i+1)-hstep)) + 2*hstep*(sin(x(i+1)+hstep) + sin(x(i+1)-hstep) - 2*sin(x(i+1)))) ...
%                  + 4*hstep*h^2*(sin(x(i+1)+hstep) + sin(x(i+1)-hstep) - 2*sin(x(i+1))) + 8*h*fi + 2*h^2*fi*(2*hstep + sin(x(i+1)+hstep) - sin(x(i+1)-hstep)));
% 
%         if i < n
%             xip2 = x(i+3);
%             fip1 = 2*xip1 + (h^2)*(xip1 + sin(xip1))-xi - xip2;
%             %g(i) = g(i) + 0.5*hstep - fip1;
%             g(i) = g(i) - fip1;
%         end
%     end
% end
