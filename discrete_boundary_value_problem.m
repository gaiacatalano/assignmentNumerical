function[F, gradf, Hess] = discrete_boundary_value_problem(n, x_bar)
    %Fun = 0;
    grad = zeros(n,1);
    H = zeros(n,n);
    h = 1/(n+1);

    x_first = 0;
    x_last = 0;
    x = [x_first; x_bar; x_last];

    f = zeros(n,1);
    for k=1:n
        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
        f(k)= (2*xk-xkm1-xkp1 + (h^2 * (xk + k*h + 1)^3)/2)^2;

        %Fun = Fun + fk^2;

        % Gradient
        df_dxkm1 = -1;
        df_dxk   = 2 + 3/2 * h^2 * (xk + k*h + 1)^2;
        df_dxkp1 = -1;

        if k > 1
            grad(k-1) = grad(k-1) + 2*(f(k) * df_dxkm1);
        end
        grad(k) = grad(k) + 2*(f(k) * df_dxk);
        if k < n
            grad(k+1) = grad(k+1) + 2*(f(k) * df_dxkp1);
        end
        
        % %Hessian
        % grad_fk = zeros(n,1);
        % if k > 1
        %     grad_fk(k-1) = df_dxkm1; 
        % end
        % grad_fk(k) = df_dxk;
        % if k < n 
        %     grad_fk(k+1) = df_dxkp1; 
        % end
        % 
        % H = H + 2*(grad_fk * grad_fk');
        % H(k,k) = H(k,k) + 2 * (f(k) * (3*h^2 * (xk + k*h + 1)));

        % Inizializzazione
        H = sparse(n, n); % Crea una matrice Hessiana sparsa n x n
        
        % Vettore gradiente parziale
        % grad_fk è un vettore colonna con al massimo 3 elementi diversi da zero
        if k > 1
            H = H + 2 * sparse([k-1; k; k+1], [k-1; k; k+1], ...
                [df_dxkm1; df_dxk; df_dxkp1] * [df_dxkm1, df_dxk, df_dxkp1], n, n);
        elseif k < n
            H = H + 2 * sparse([k; k+1], [k; k+1], ...
                [df_dxk; df_dxkp1] * [df_dxk, df_dxkp1], n, n);
        else
            H = H + 2 * sparse(k, k, df_dxk^2, n, n);
        end
        % Aggiunta del termine diagonale specifico
        H(k,k) = H(k,k) + 2 * (f(k) * (3*h^2 * (xk + k*h + 1)));
    end

    H = 0.5 * (H + H');  % forza simmetria numerica
    
    F = @(x) sum(f.^2);

    gradf = @(x) grad;

    Hess = @(x) H;
end