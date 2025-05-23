function H = problem_213_hess_fd(x, hstep)
% per l'hessiana uso la differenza centrata prchè più precisa per stimare
% la curvatura

    if nargin < 2
        hstep = 1e-5;
    end

    n = length(x);
    H = sparse(n,n);  % la Hessiana è sparsa
    h = 1 / (n + 1);

    d0 = zeros(n,1);
    d1 = zeros(n,1);
    d2 = zeros(n,1);

    % Estendi x con condizioni al contorno
    x_ext = [0; x; 1];

    % Calcola il gradiente originale (strutturato)
    %g0 = problem_213_grad_df(x, hstep);  % usa il gradiente strutturato già implementato

    % Ciclo sulle colonne (variabili x_j)
    for i = 1:n
        xim1 = x_ext(i);
        xi = x_ext(i+1);
        xip1 = x_ext(i+2);

        if i > 1
            d0(i) = d0(i) + 2;
        end

        d0(i) = d0(i) + 8 + 8*(h^2) + 2*(h^4)*(hstep^2 + sin(xi)^2 + 2*cos(xi)*sin(hstep)/hstep) -(4/hstep^2)*(sin(xi)^2)*cos(hstep) + (8*(h^2)/hstep)*cos(xi)*sin(hstep); 

        if i < n
            d0(i) = d0(i) + 2;
            d1(i) = d1(i) +  (-2*(h^2)/hstep)*(hstep +sin(xi + hstep)-sin(xi)) -4 - (-2*(h^2)/hstep)*(hstep+sin(xip1 + hstep)-sin(xip1));    %forse sommarlo non serve, lo riempio una volta sola
        end

        if i < n-1
            d2(i) = 2;
        end

    end

    H = spdiags([d2 d1 d0 d1 d2], [-2 -1 0 -1 -2], n, n);

    % Simmetrizza per sicurezza
    H = 0.5 * (H + H');
end