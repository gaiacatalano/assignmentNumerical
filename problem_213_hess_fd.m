function H = problem_213_hess_fd(x, hstep)
% per l'hessiana uso la differenza centrata prchè più precisa per stimare
% la curvatura

    if nargin < 2
        hstep = 1e-5;
    end

    n = length(x);
    H = sparse(n,n);  % la Hessiana è sparsa
    h = 1 / (n + 1);

    % Estendi x con condizioni al contorno
    %x_ext = [0; x; 1];

    % Calcola il gradiente originale (strutturato)
    %g0 = problem_213_grad_df(x, hstep);  % usa il gradiente strutturato già implementato

    % Ciclo sulle colonne (variabili x_j)
    for i = 1:n
        % Perturbazione positiva
        xp = x;
        xp(i) = xp(i) + hstep;
        gp = problem_213_grad_fd(xp, hstep);

        % Perturbazione negativa
        xm = x;
        xm(i) = xm(i) - hstep;
        gm = problem_213_grad_fd(xm, hstep);

        % Derivata centrale ∂grad_i/∂x_j = colonna j della Hessiana
        H(:,i) = (gp - gm) / (2*hstep);
    end

    % Simmetrizza per sicurezza
    H = 0.5 * (H + H');
end