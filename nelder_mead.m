function[simplex, k] = nelder_mead(x0, fun, rho, chi, gamma, sigma, kmax, tol)

    % Parametro di perturbazione
    delta = 0.05 * norm(x0);

    % Più è grande delta scelto, più iterazioni sono necessarie alla
    % convergenza

    % Genero gli altri punti del simplesso
    x1 = x0 + [delta; 0];
    x2 = x0 + [0; delta];

    % Creo il simplesso
    simplex = [x0, x1, x2];

    f = zeros(1,3);
  
    % Valuto la funzione nei punti del simplesso
    for i = 1:3
        f(i) = fun(simplex(:,i));
    end
    prev_f = f;

    for k=1:kmax
        
        % Ordino i valori crescenti e ottengo anche gli indici
        [f_sorted, idx] = sort(f);

        % Riordino il simplesso in base agli indici
        simplex = simplex(:, idx);

        % Stampa i punti del simplesso ordinato e i valori della funzione
        disp('Punti del simplesso ordinato e valori f:');
        for i = 1:3
            fprintf('x%d = [', i);
            fprintf('%f ', simplex(:,i));
            fprintf('], f = %.6f\n', f_sorted(i));
        end

        
        % Baricentro di n punti migliori
        x_bar = mean(simplex(:,1:2), 2);
        
        % Riflessione
        x_r = x_bar + rho*(x_bar - simplex(:,end));
        f_r = fun(x_r);

        % se f_r compresa tra f_1 e f_n -> accetto x_r al posto di x_n+1
        if(f_sorted(1) <= f_r && f_r <= f_sorted(2))
           simplex(:,3) = x_r;
           disp("Accetto riflessione")
 
        % altrimenti se f_r < f_1 -> espansione
        elseif(f_r < f_sorted(1))
            x_e = x_bar + chi*(x_r - x_bar);
            f_e = fun(x_e);

            if(f_e < f_r)
                simplex(:,3) = x_e;
                %continue
                disp("Accetto espansione")
            
            else
                simplex(:,3) = x_r;
                %continue
                disp("Accetto riflessione non espansione")

            end
    
        else % contrazione
            x_c = x_bar - gamma*(x_bar - x_r);
            f_c = fun(x_c);
            if(f_c < f_sorted(3))
                simplex(:,3) = x_c;
                disp("Accetto contrazione")

            else % shrink
                simplex(:,2) = simplex(:,1) + sigma * (simplex(:,2) - simplex(:,1));
                simplex(:,3) = simplex(:,1) + sigma * (simplex(:,3) - simplex(:,1));
                disp("Accetto shrink")
    
            end

        end


    % Condizioni di arresto basate su diametro e differenze dei valori di f
    % calcolata nel simplesso attuale

    for i = 1:3
            f(i) = fun(simplex(:,i));
    end
    
    fprintf('Iter %d: x = [%f, %f], f = %.6f\n', k, simplex(1,1), simplex(2,1), f(1));
    fprintf('Iter %d: x = [%f, %f], f = %.6f\n', k, simplex(1,2), simplex(2,2), f(2));
    fprintf('Iter %d: x = [%f, %f], f = %.6f\n', k, simplex(1,3), simplex(2,3), f(3));
    
    % Criterio di arresto basato sulla variazione dei valori della
    % funzione, diametro del simplesso, valore massimo e minimo della
    % funzione

    if max(abs(f - prev_f)) < tol && norm(simplex(2:end, :)-simplex(1,:)) < tol &&  abs(max(f) - min(f)) < tol
        disp('Criterio di arresto attivato per Nelder-Mead');
        break;
    end

    prev_f = f;  % Memorizza i valori per la prossima iterazione

    end

end



