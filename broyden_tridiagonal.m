function[F, gradf, Hess] = broyden_tridiagonal(n, x_bar)
    %F=0;
    grad = zeros(n,1);
    H = zeros(n,n);
    f = zeros(n,1);

    x_first = 0;
    x_last = 0;
    x = x_bar*ones(n,1);
    x = [x_first; x; x_last];
   
    for k=1:n
        idx = [];
        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);

        f(k) = (3-2*xk)*xk - xkm1 -2*xkp1 + 1;
        %F  = F + 0.5 * fk^2;
        
        % Gradient
        df_dxkm1 = -1;
        df_dxk   = 3 - 4*xk;
        df_dxkp1 = -2;

        if k > 1
            grad(k-1) = grad(k-1) + f(k) * df_dxkm1;
            idx(end+1) = k - 1;
        end
        grad(k) = grad(k) + f(k) * df_dxk;
        idx(end+1) = k;
        if k < n
            grad(k+1) = grad(k+1) + f(k) * df_dxkp1;
            idx(end+1) = k + 1;
        end
        
        %Hessian
        % Inizializzazione Hessiana sparsa
        H = sparse(n, n);  % matrice n x n
        grad_fk = zeros(n,1);
        if k > 1
            grad_fk(k-1) = df_dxkm1; 
        end
        grad_fk(k) = df_dxk;
        if k < n 
            grad_fk(k+1) = df_dxkp1; 
        end
        
        % Nota: questo aggiorna solo gli indici idx
        H = H + sparse(idx, idx', grad_fk * grad_fk', n, n);
        
        % Aggiunta del termine diagonale costante
        H(k, k) = H(k, k) + f(k) * (-4);  % d²f_k/dx_k² = -4
        
        % if k > 1
        %     grad_fk(k-1) = df_dxkm1; 
        % end
        % grad_fk(k) = df_dxk;
        % if k < n 
        %     grad_fk(k+1) = df_dxkp1; 
        % end
        % 
        % H = H + grad_fk * grad_fk';
        % H(k,k) = H(k,k) + f(k) * (-4);  % d²f_k/dx_k² = -4
        % 

        % versione sparsa dell'hessiana
        % rows = []; cols = []; vals=[];
        % if k > 1
        %     rows = [rows; k-1; k-1; k-1];
        %     cols = [cols; k-1; k;   k+1];
        %     vals = [vals; df_dxkm1; df_dxkm1*df_dxk; df_dxkm1*df_dxkp1];
        % end
        % 
        % rows = [rows; k; k; k];
        % cols = [cols; k; k-1*(k>1); k+1*(k<n)];
        % vals = [vals; df_dxk^2 + f(k)*(-4); ...
        %                df_dxk*df_dxkm1*(k>1); ...
        %                df_dxk*df_dxkp1*(k<n)];
        % 
        % if k < n
        %     rows = [rows; k+1; k+1; k+1];
        %     cols = [cols; k+1; k-1; k];
        %     vals = [vals; df_dxkp1^2; df_dxkp1*df_dxkm1*(k>1); df_dxkp1*df_dxk];
        % end
    end

    % Costruzione Hessiana sparsa simmetrica
    % H = sparse(rows, cols, vals, n, n);
    % H = 0.5 * (H + H');  % forza simmetria numerica

    F = @(x) 0.5 * sum(f.^2);

    gradf = @(x) grad;

    Hess = @(x) H;
end