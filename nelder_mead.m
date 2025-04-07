function[simplex] = nelder_mead(x0, fun, n, rho, chi, gamma, sigma, kmax, tol_x, max_no_improvement)

    % Parametro di perturbazione
    delta = 0.1;

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

    for k=1:kmax
        
        % Ordino i valori crescenti e ottengo anche gli indici
        [f_sorted, idx] = sort(f);
        
        % Riordino il simplesso in base agli indici
        simplex = simplex(:, idx);
        
        % Baricentro di n punti migliori
        x_bar = mean(simplex(:,1:n), 2);
        
        % Riflessione
        x_r = x_bar + rho*(x_bar - x(idx(n+1)));
        f_r = f(x_r);

        % se f_r compresa tra f_1 e f_n -> accetto x_r al posto di x_n+1
        if(f_sorted(1) <= f_r && f_r > f_sorted(n))
           simplex(:,n+1) = x_r;
           continue
 
        % altrimenti se f_r < f_1 -> espansione
        elseif(f_r < f_sorted(1))
            x_e = x_bar + chi*(x_r - x_bar);
            f_e = f(x_e);

            if(f_e < f_r)
                simplex(:,n+1) = x_e;
                continue
            
            else
                simplex(:,n+1) = x_r;
                continue
            end
    
        else % contrazione
            x_c = x_bar - gamma*(x_bar - x_r);
            f_c = f(x_c);
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

    diametro = max(vecnorm(simplex - mean(simplex, 2), 2, 1));
    if diametro < tol_x || (max(f) - min(f) < tol_f)
        no_improvement_count = no_improvement_count + 1;
    else
        no_improvement_count = 0;  % Reset se c'è stato miglioramento
    end

    % Verifica il numero di iterazioni senza miglioramenti
    if no_improvement_count >= max_no_improvement
        disp('Nessun miglioramento significativo, fermo!');
        break;
    end

    end

end



