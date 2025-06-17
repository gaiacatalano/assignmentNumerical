function[simplex] = nelder_mead_n(x0, fun, n, rho, chi, gamma, sigma, kmax, tol)

    % Perturbation parameter
    delta = 0.12 * norm(x0);

    % Generation of simplex points
    simplex = x0.*ones(1, n+1);
    for k=1:n
        simplex(k,k+1) = simplex(k,k+1) + delta;
    end

    f = zeros(1,n+1);

    for i = 1:n+1
        f(i) = fun(simplex(:,i));
    end
    prev_f = f;

    for k=1:kmax
        
        % Sort the values in ascending order and get the indices
        [f_sorted, idx] = sort(f);

        % Reorder the simplex
        simplex = simplex(:, idx);

        % Centroid of n best points
        x_bar = mean(simplex(:,1:n), 2);
        
        % Reflection
        x_r = x_bar + rho*(x_bar - simplex(:,end));
        f_r = fun(x_r);

        % If f_r in (f_1,f_n) -> accept x_r
        if(f_sorted(1) <= f_r && f_r <= f_sorted(n))
           simplex(:,n+1) = x_r;
 
        % Otherwise -> Expansion
        elseif(f_r < f_sorted(1))
            x_e = x_bar + chi*(x_r - x_bar);
            f_e = fun(x_e);

            if(f_e < f_r)
                simplex(:,n+1) = x_e;
            else
                % Accept reflection
                simplex(:,n+1) = x_r;
            end
    
        else % Contraction
            x_c = x_bar - gamma*(x_bar - x_r);
            f_c = fun(x_c);
            if(f_c < f_sorted(n+1))
                simplex(:,n+1) = x_c;

            else % Shrink
                for i=2:n+1
                    simplex(:,i) = simplex(:,1) + sigma * (simplex(:,i)-simplex(:,1));
                end
    
            end

        end

    for i = 1:n+1
            f(i) = fun(simplex(:,i));
    end
    
    % Stopping criterion based on function value variation, 
    % simplex diameter, and the function's maximum and minimum values

    if max(abs(f - prev_f)) < tol && norm(simplex(2:end, :)-simplex(1,:)) < tol &&  abs(max(f) - min(f)) < tol
        break;
    end

    prev_f = f;  % Save for next iteration

    end

end



