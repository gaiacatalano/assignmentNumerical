function H = problem_213_hess_fd(x, hstep, bool_hstep_i)
% per l'hessiana uso la differenza centrata prchè più precisa per stimare
% la curvatura

    if nargin < 2
        hstep = 1e-5;
    end

    if bool_hstep_i==1
        hstep_i= abs(x)*hstep;
    end

    n = length(x);
    %H = sparse(n,n);  % la Hessiana è sparsa
    h = 1 / (n + 1);

    d0 = zeros(n,1);
    %dp1 = zeros(n,1);
    dm1 = zeros(n,1);
    %dp2 = zeros(n,1);
    dm2 = zeros(n,1);

    % Estendi x con condizioni al contorno
    x = [0; x; 1];

    % Calcola il gradiente originale (strutturato)
    %g0 = problem_213_grad_df(x, hstep);  % usa il gradiente strutturato già implementato

    % Ciclo sulle colonne (variabili x_j)
    for i = 1:n

        fi = 2*x(i+1) + h^2*(x(i+1) + sin(x(i+1))) - x(i) - x(i+2);

        if bool_hstep_i == 0
            d0(i) = 4 + (h^4)/(2*hstep^2)*(2*hstep^2 + sin(x(i+1) + hstep)^2 + sin(x(i+1)-hstep)^2 + 2*sin(x(i+1))^2 - 2*sin(x(i+1)+hstep)*sin(x(i+1)) - 2*sin(x(i+1)-hstep)*sin(x(i+1)) + 2*hstep*(sin(x(i+1)+hstep) - sin(x(i+1)-hstep))) ...
                + (h^2)*(2/hstep*(2*hstep + sin(x(i+1)+hstep) - sin(x(i+1)-hstep)) + fi/(hstep^2) * (sin(x(i+1)+hstep) + sin(x(i+1)-hstep)-2*sin(x(i+1))));
        else
            d0(i) = 4 + (h^4) / (2 * hstep_i(i)^2) * (2 * hstep_i(i)^2 + sin(x(i+1) + hstep_i(i))^2 + sin(x(i+1) - hstep_i(i))^2 + 2 * sin(x(i+1))^2 - 2 * sin(x(i+1) + hstep_i(i)) * sin(x(i+1)) - 2 * sin(x(i+1) - hstep_i(i)) * sin(x(i+1)) + 2 * hstep_i(i) * (sin(x(i+1) + hstep_i(i)) - sin(x(i+1) - hstep_i(i)))) ...
                + h^2 * ((2 / hstep_i(i)) * (2 * hstep_i(i) + sin(x(i+1) + hstep_i(i)) - sin(x(i+1) - hstep_i(i))) + (fi / hstep_i(i)^2) * (sin(x(i+1) + hstep_i(i)) + sin(x(i+1) - hstep_i(i)) - 2 * sin(x(i+1))));
        end

        %d0(i) = 4 + h^4*(hstep^2 + sin(x(i+1))^2*cos(hstep)^2 + cos(x(i+1))^2*sin(hstep)^2 + sin(x(i+1))^2 +2*hstep*cos(x(i+1))*sin(hstep) ...
        %            - 2*sin(x(i+1))^2*cos(hstep))/(hstep^2) + 4*h^2*(hstep + cos(x(i+1))*sin(hstep))/hstep + 2*fi*(h^2*(sin(x(i+1))*cos(h) - sin(x(i+1))))/(hstep^2);
        %plus = h^4/2*((hstep + sin(x(i+1) + hstep) - sin(x(i+1)))^2)/(hstep^2) - h^2*(hstep + sin(x(i+1) + hstep) - sin(x(i+1)))/hstep ...
        %            + h^2*(fi*(hstep + sin(x(i+1) + hstep) - sin(x(i+1))))/(hstep^2);
        
        if i ~= 1
            d0(i) = d0(i) + 1;

            % if bool_hstep_i==0
            %     dm1(i-1) = dm1(i-1) + plus;
            % else
            %     dm1(i-1)
            % end
        end
            
        if i ~= n
            d0(i) = d0(i) + 1;

            if bool_hstep_i==0
                dm1(i) = (h^2)/(hstep)*(sin(x(i+1)) + sin(x(i+2)) - sin(x(i+1)+hstep) - sin(x(i+2)+hstep) - 2*hstep) -4;
            else
                dm1(i) = (h^2)/(hstep_i(i))*(sin(x(i+1)) - sin(x(i+1)+hstep_i(i)) - hstep_i(i)) + (h^2)/(hstep_i(i+1))*(sin(x(i+2)) - sin(x(i+2)+hstep_i(i+1)) - hstep_i(i+1)) - 4;
            end
        end
    
        if i < n-1
            dm2(i) = 1;
        end

        % xim1 = x_ext(i);
        % xi = x_ext(i+1);
        % xip1 = x_ext(i+2);
        % 
        % if i > 1
        %     d0(i) = d0(i) + 2;
        % 
        % end
        % 
        % d0(i) = d0(i) + 8 + 8*(h^2) + 2*(h^4)*(hstep^2 + sin(xi)^2 + 2*cos(xi)*sin(hstep)/hstep) -(4/hstep^2)*(sin(xi)^2)*cos(hstep) + (8*(h^2)/hstep)*cos(xi)*sin(hstep); 
        % 
        % if i < n
        %     d0(i) = d0(i) + 2;
        %     dm1(i) = dm1(i) +  (-2*(h^2)/hstep)*(hstep +sin(xi + hstep)-sin(xi)) -4 - (-2*(h^2)/hstep)*(hstep+sin(xip1 + hstep)-sin(xip1));    %forse sommarlo non serve, lo riempio una volta sola
        %     %dp1(i+1) = dp1(i+1) +  (-2*(h^2)/hstep)*(hstep +sin(xi + hstep)-sin(xi)) -4 - (-2*(h^2)/hstep)*(hstep+sin(xip1 + hstep)-sin(xip1));    %forse sommarlo non serve, lo riempio una volta sola
        % end
        % 
        % if i < n-1
        %     dm2(i) = 2;
        %     %dp2(i+2) = 2;
        % 
        % end

    end

    H = spdiags([dm2 dm1 d0 [0;dm1(1:end-1)] [0;0;dm2(1:end-2)]], [-2 -1 0 -1 -2], n, n);

    % Simmetrizza per sicurezza
    H = 0.5 * (H + H');
    %any(isnan(H(:))) || any(isinf(H(:)))

end