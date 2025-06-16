function H = discrete_boundary_value_hess_fd(x, hstep, bool_hstep_i)
    % Approximate Hessian for the discrete BVP via finite differences
    %
    % Input:
    %   x     - n x 1 vector
    %   hstep - finite difference step (default: 1e-5)
    %
    % Output:
    %   H     - n x n sparse Hessian approximation

    if nargin < 2
        hstep = 1e-5;
    end

    if bool_hstep_i==1
        hstep_i= abs(x)*hstep;
    end

    n = length(x);
    %H = sparse(n, n);
    x_ext = [0; x; 0];
    h = 1 / (n + 1);
    d0 = zeros(n,1);
    %dp1 = zeros(n,1);
    dm1 = zeros(n,1);
    %dp2 = zeros(n,1);
    dm2 = zeros(n,1);

    for i = 1:n
        xim1 = x_ext(i);
        xi = x_ext(i+1);
        xip1 = x_ext(i+2);

        if bool_hstep_i==1
            hstep= hstep_i(i);
        end

        fi = 2*xi - xim1 - xip1 + (h^2 / 2)*(xi + i*h + 1)^3;

        d0(i) = d0(i) + 2*(2+((h*hstep)^2)/2)^2 + 6*fi*(h^2)*(xi+i*h + 1) + ((9*h^4)/2)*((xi + i*h + 1)^2)*((xi+i*h + 1)^2 + hstep^2) + 3*(h^2)*(4+(h*hstep)^2)*(xi + i*h+1)^2;

        if i ~= 1
            d0(i) = d0(i) + 2;
            
        end 

        if i ~= n
            %xip2 = x_ext(i+3);
            %fip1 = 2*xip1 - xi - xip2 + (h^2 / 2)*(xip1 + (i+1)*h + 1)^3;
            d0(i) = d0(i) + 2;
            if bool_hstep_i == 0
                dm1(i) = -8-3*(h^2)*(xi+i*h+1)*((xi+i*h+1)+hstep)-(h*hstep)^2 -(h*hstep)^2 - 3*(h^2)*(xip1+(i+1)*h+1)*((xip1+(i+1)*h+1)+hstep);
                %dp1(i+1) = dp1(i+1) -3*(h^2)*(xi+i*h+1)*((xi+i*h+1)+hstep)-(h*hstep)^2 -4 -(h*hstep)^2 - 3*(h^2)*(xip1+(i+1)*h+1)*((xip1+(i+1)*h+1)+hstep);
            else
                dm1(i) = -2*(hstep_i(i+1) * (2*hstep_i(i) + (h^2/2)*(hstep_i(i)^3 + 3*hstep_i(i)*(xi + i*h + 1)*(hstep_i(i) + xi + i*h +1)))  ...
                    + hstep_i(i)*(2*hstep_i(i+1) + (h^2/2)*(hstep_i(i+1)^3 + 3*hstep_i(i+1)*(xip1 + (i+1)*h + 1)*(hstep_i(i+1) + xip1 + (i+1)*h +1))));
            end
        end 

        if i < n-1
            dm2(i) = 2;
            %dp2(i+2) = 2;
        end

    end

    H = spdiags([dm2 dm1 d0 [0; dm1(1:end-1)] [0;0;dm2(1:end-2)]], [-2 -1 0 -1 -2], n, n);

    % Enforce symmetry
    H = 0.5 * (H + H');
end