% Approximates the Hessian of the Discrete Boundary Value problem 
% using finite difference methods

function H = discrete_boundary_value_hess_fd(x, hstep, bool_hstep_i)

    if nargin < 2 || isempty(hstep)
        hstep = 1e-5;
    end

    if nargin < 3 || isempty(bool_hstep_i)
        bool_hstep_i = 0;
    end

    n = length(x);
    H = sparse(n,n);

    if bool_hstep_i==1
        hstep_i= abs(x)*hstep;
    else
        hstep_i = hstep * ones(n, 1);
    end

    %x_ext = [0; x; 0];
    %h = 1 / (n + 1);

    % d0 = zeros(n,1);
    % dm1 = zeros(n,1);
    % dm2 = zeros(n,1);

    for i = 1:n

        hi = hstep_i(i);
        e = zeros(n,1);
        e(i) = 1;

        % gradiente analitico nei due punti perturbati
        gp  = discrete_boundary_value_grad(x + hi*e);
        gm = discrete_boundary_value_grad(x - hi*e);

        % i-esima colonna della Hessiana (d/dx_i del gradiente)
        H(:, i) = (gp - gm) / (2 * hi);
    end
    
    H = 0.5 * (H + H.');
    
end


%         xim1 = x_ext(i);
%         xi = x_ext(i+1);
%         xip1 = x_ext(i+2);
% 
%         if bool_hstep_i==1
%             hstep= hstep_i(i);
%         end
% 
%         fi = 2*xi - xim1 - xip1 + (h^2 / 2)*(xi + i*h + 1)^3;
% 
%         d0(i) = d0(i) + 2*(2+((h*hstep)^2)/2)^2 + 6*fi*(h^2)*(xi+i*h + 1) + ...
%         ((9*h^4)/2)*((xi + i*h + 1)^2)*((xi+i*h + 1)^2 + hstep^2) + ... 
%         3*(h^2)*(4+(h*hstep)^2)*(xi + i*h+1)^2;
% 
%         if i ~= 1
%             d0(i) = d0(i) + 2;
% 
%         end 
% 
%         if i ~= n
%             d0(i) = d0(i) + 2;
% 
%             if bool_hstep_i == 0
%                 dm1(i) = -8-3*(h^2)*(xi+i*h+1)*((xi+i*h+1)+hstep)-(h*hstep)^2 ...
%                 -(h*hstep)^2 - 3*(h^2)*(xip1+(i+1)*h+1)*((xip1+(i+1)*h+1)+hstep);
%             else
%                 dm1(i) = -2*(hstep_i(i+1) * (2*hstep_i(i) + (h^2/2)*(hstep_i(i)^3 + ...
%                     3*hstep_i(i)*(xi + i*h + 1)*(hstep_i(i) + xi + i*h +1)))  ...
%                     + hstep_i(i)*(2*hstep_i(i+1) + (h^2/2)*(hstep_i(i+1)^3 + ...
%                     3*hstep_i(i+1)*(xip1 + (i+1)*h + 1)*(hstep_i(i+1) + xip1 + (i+1)*h +1))));
%             end
% 
%         end 
% 
%         if i < n-1
%             dm2(i) = 2;
%         end
% 
%     end
% 
%     H = spdiags([dm2 dm1 d0 [0; dm1(1:end-1)] [0;0;dm2(1:end-2)]], [-2 -1 0 -1 -2], n, n);
%     H = 0.5 * (H + H');
% 
% end