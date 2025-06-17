% Approximates the Hessian of the Problem 213
% using finite difference methods

function H = problem_213_hess_fd(x, hstep, bool_hstep_i)

    if nargin < 2
        hstep = 1e-5;
    end

    if bool_hstep_i==1
        hstep_i= abs(x)*hstep;
    end

    n = length(x);
    h = 1 / (n + 1);

    d0 = zeros(n,1);
    dm1 = zeros(n,1);
    dm2 = zeros(n,1);

    x = [0; x; 1];

    for i = 1:n

        fi = 2*x(i+1) + h^2*(x(i+1) + sin(x(i+1))) - x(i) - x(i+2);

        if bool_hstep_i == 0
            d0(i) = 4 + (h^4)/(2*hstep^2)*(2*hstep^2 + sin(x(i+1) + hstep)^2 + ...
                sin(x(i+1)-hstep)^2 + 2*sin(x(i+1))^2 - 2*sin(x(i+1)+hstep)*sin(x(i+1)) -...
                2*sin(x(i+1)-hstep)*sin(x(i+1)) + 2*hstep*(sin(x(i+1)+hstep) - ...
                sin(x(i+1)-hstep))) + (h^2)*(2/hstep*(2*hstep + sin(x(i+1)+hstep) - ...
                sin(x(i+1)-hstep)) + fi/(hstep^2) * (sin(x(i+1)+hstep) + sin(x(i+1)-hstep)-...
                2*sin(x(i+1))));
        else
            d0(i) = 4 + (h^4) / (2 * hstep_i(i)^2) * (2 * hstep_i(i)^2 + sin(x(i+1) +...
                hstep_i(i))^2 + sin(x(i+1) - hstep_i(i))^2 + 2 * sin(x(i+1))^2 - ...
                2 * sin(x(i+1) + hstep_i(i)) * sin(x(i+1)) - 2 * sin(x(i+1) - hstep_i(i)) * sin(x(i+1)) + ...
                2 * hstep_i(i) * (sin(x(i+1) + hstep_i(i)) - sin(x(i+1) - hstep_i(i)))) ...
                + h^2 * ((2 / hstep_i(i)) * (2 * hstep_i(i) + sin(x(i+1) + hstep_i(i)) ...
                - sin(x(i+1) - hstep_i(i))) + (fi / hstep_i(i)^2) * (sin(x(i+1) + hstep_i(i)) + sin(x(i+1) - hstep_i(i)) ...
                - 2 * sin(x(i+1))));
        end

        if i ~= 1
            d0(i) = d0(i) + 1;
        end
            
        if i ~= n
            d0(i) = d0(i) + 1;

            if bool_hstep_i==0
                dm1(i) = (h^2)/(hstep)*(sin(x(i+1)) + sin(x(i+2)) - sin(x(i+1)+hstep) -...
                    sin(x(i+2)+hstep) - 2*hstep) -4;
            else
                dm1(i) = (h^2)/(hstep_i(i))*(sin(x(i+1)) - sin(x(i+1)+hstep_i(i)) - hstep_i(i)) +...
                    (h^2)/(hstep_i(i+1))*(sin(x(i+2)) - sin(x(i+2)+hstep_i(i+1)) - hstep_i(i+1)) - 4;
            end
        end
    
        if i < n-1
            dm2(i) = 1;
        end   

    end

    H = spdiags([dm2 dm1 d0 [0;dm1(1:end-1)] [0;0;dm2(1:end-2)]], [-2 -1 0 -1 -2], n, n);
    H = 0.5 * (H + H');

end