function g = problem_213_grad_fd(x, hstep)
    if nargin < 2
        hstep = 1e-6;
    end

    n = length(x);
    h = 1 / (n + 1);  
    g = zeros(n,1);
    x_ext = [0; x; 0];  % estensione con x0 = 0, x_{n+1} = 0

    for i = 1:n
        fx=0;

        if i > 1
            xim1 = x_ext(i-1);
            xi   = x_ext(i);
            xip1 = x_ext(i+1);
            fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
            fx = fx + fi^2;            
        end
        
        xim1 = x_ext(i);
        xi   = x_ext(i+1);
        xip1 = x_ext(i+2);
        fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
        fx = fx + fi^2;

        if i < n
            xim1 = x_ext(i+1);
            xi   = x_ext(i+2);
            xip1 = x_ext(i+3);
            fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
            fx = fx + fi^2;
        end
        
        %pezzo perturbato
        xh = x;
        xh(i) = xh(i) + hstep;
        xh_ext = [0; xh; 0];

        fxh = 0;
        if i > 1
            xim1 = xh_ext(i-1);
            xi   = xh_ext(i);
            xip1 = xh_ext(i+1);
            fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
            fxh = fxh + fi^2;            
        end
        
        xim1 = xh_ext(i);
        xi   = xh_ext(i+1);
        xip1 = xh_ext(i+2);
        fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
        fxh = fxh + fi^2;

        if i < n
            xim1 = xh_ext(i+1);
            xi   = xh_ext(i+2);
            xip1 = xh_ext(i+3);
            fi = 2*xi + (h^2)*(xi + sin(xi)) - xim1 - xip1;
            fxh = fxh + fi^2;
        end

        g(i) = (fxh - fx) / hstep;

    end


end
