function[simplex] = nelder_mead(x0, fun, n, rho, chi, gamma, sigma, kmax, tol_x, max_no_improvement, tol_f)

    % Parametro di perturbazione
    delta = 0.05;

    no_improvement_count = 0;

    % Genero gli altri punti del simplesso
    % GENERALIZZARE CON N SE POSSIBILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    x1 = x0 + [delta; 0];
    x2 = x0 + [0; delta];

    % Creo il simplesso
    simplex = [x0, x1, x2];

    f = zeros(1,n+1);
  
    % Valuto la funzione nei punti del simplesso
    for i = 1:n+1
            f(i) = fun(simplex(:,i));
    end
    prev_f = f;

    for k=1:kmax
        
        % Ordino i valori crescenti e ottengo anche gli indici
        [f_sorted, idx] = sort(f);
        
        % Riordino il simplesso in base agli indici
        simplex = simplex(:, idx);
        
        % Baricentro di n punti migliori
        x_bar = mean(simplex(:,1:n), 2);
        
        % Riflessione
        x_r = x_bar + rho*(x_bar - simplex(:,end));
        f_r = fun(x_r);

        % se f_r compresa tra f_1 e f_n -> accetto x_r al posto di x_n+1
        if(f_sorted(1) <= f_r && f_r > f_sorted(n))
           simplex(:,n+1) = x_r;
           %continue
 
        % altrimenti se f_r < f_1 -> espansione
        elseif(f_r < f_sorted(1))
            x_e = x_bar + chi*(x_r - x_bar);
            f_e = fun(x_e);

            if(f_e < f_r)
                simplex(:,n+1) = x_e;
                %continue
            
            else
                simplex(:,n+1) = x_r;
                %continue
            end
    
        else % contrazione
            x_c = x_bar - gamma*(x_bar - x_r);
            f_c = fun(x_c);
            if(f_c < f_sorted(n+1))
                simplex(:,n+1) = x_c;
            else % shrink % SISTEMARE PER GENERALIZZARE !!!!!!!
                simplex(:,2) = simplex(:,1) + sigma * (simplex(:,2) - simplex(:,1));
                simplex(:,3) = simplex(:,1) + sigma * (simplex(:,3) - simplex(:,1));
            end

        end


    % Condizioni di arresto basate su diametro e differenze dei valori di f
    % calcolata nel simplesso attuale

    for i = 1:n+1
            f(i) = fun(simplex(:,i));
    end

    fprintf('Iter %d: x = [%f, %f], f = %.6f\n', k, simplex(1,1), simplex(2,1), f(1));

    % Aggiorna il contatore dei miglioramenti
    if max(abs(f - prev_f)) < tol_f
        no_improvement_count = no_improvement_count + 1;
    else
        no_improvement_count = 0;
    end

    % Criterio di arresto basato sulla variazione dei valori della funzione
    if max(abs(f - prev_f)) < tol_f || max(vecnorm(simplex - mean(simplex, 2), 2, 1)) < tol_x
        disp('Criterio di arresto attivato per Nelder-Mead');
        break;
    end

    
    prev_f = f;  % Memorizza i valori per la prossima iterazione

    % diametro = max(vecnorm(simplex - mean(simplex, 2), 2, 1));
    % if diametro < tol_x || (max(f) - min(f) < tol_f)
    %     no_improvement_count = no_improvement_count + 1;
    % else
    %     no_improvement_count = 0;  % Reset se c'è stato miglioramento
    % end

    % Verifica il numero di iterazioni senza miglioramenti
    if no_improvement_count >= max_no_improvement
        disp('Nessun miglioramento significativo, fermo!');
        break;
    end
    %fprintf('Iter %d: Best f = %.6f at x = [%f, %f]\n', k, f_sorted(1), simplex(1,1), simplex(2,1));

    end

end



