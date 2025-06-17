function[simplex, k] = nelder_mead(x0, fun, rho, chi, gamma, sigma, kmax, tol)

    % Perturbation parameter
    delta = 0.05 * norm(x0);

    % Generation of simplex points
    x1 = x0 + [delta; 0];
    x2 = x0 + [0; delta];

    simplex = [x0, x1, x2];

    f = zeros(1,3);
  
    for i = 1:3
        f(i) = fun(simplex(:,i));
    end
    prev_f = f;

    for k=1:kmax
        
        % Sort the values in ascending order and get the indices
        [f_sorted, idx] = sort(f);

        % Reorder the simplex
        simplex = simplex(:, idx);

        % Centroid of n best points
        x_bar = mean(simplex(:,1:2), 2);
        
        % Reflection
        x_r = x_bar + rho*(x_bar - simplex(:,end));
        f_r = fun(x_r);

        % If f_r in (f_1,f_n) -> accept x_r
        if(f_sorted(1) <= f_r && f_r <= f_sorted(2))
           simplex(:,3) = x_r;
 
        % Otherwise -> Expansion
        elseif(f_r < f_sorted(1))
            x_e = x_bar + chi*(x_r - x_bar);
            f_e = fun(x_e);

            if(f_e < f_r)
                simplex(:,3) = x_e;            
            else
                % Accept reflection
                simplex(:,3) = x_r;
            end
    
        else % Contraction
            x_c = x_bar - gamma*(x_bar - x_r);
            f_c = fun(x_c);
            if(f_c < f_sorted(3))
                simplex(:,3) = x_c;

            else % Shrink
                simplex(:,2) = simplex(:,1) + sigma * (simplex(:,2) - simplex(:,1));
                simplex(:,3) = simplex(:,1) + sigma * (simplex(:,3) - simplex(:,1));    
            end

        end

    for i = 1:3
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



