function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    modified_newton_backtracking_preconditioning(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax, hstep, hstep_i)

% Armijo condition
farmijo = @(fk, alpha, c1_gradfk_pk) fk + alpha * c1_gradfk_pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
xk = x0;
fk = f(xk);

% To split the cases using finite differences
if nargin == 9
    gradfk = gradf(xk);
else
    gradfk = gradf(xk, hstep, hstep_i);
end

k = 0;
gradfk_norm = norm(gradfk);

while k < kmax && gradfk_norm >= tolgrad
    if nargin == 9
        Hessfk = Hessf(xk);
    else
        Hessfk = Hessf(xk, hstep, hstep_i);
    end
    
    beta = norm(Hessfk, 'fro');

    if min(diag(Hessfk)) > 0
        tau_k = 0;
    else
        tau_k = beta / 2;
    end

    % Regularization process
    j = 0;
    while true
        try
            R = chol(Hessfk + tau_k * speye(size(Hessfk)));
            break;
        catch
            tau_k = max(2 * tau_k, beta / 2);
            j = j + 1;
        end
    end

    E_k = tau_k*speye(size(Hessfk,1));
    B_k = Hessfk + E_k;

    % Solve with pcg with diagonal preconditioning
    M_diag = diag(B_k);
    M_diag(M_diag <= 0) = 1;
    n = length(xk);
    M = spdiags(M_diag, 0, n, n);

    [pk, ~] = pcg(B_k, -gradfk, 1e-6, kmax, M);

    alpha = 0.5;
    
    xnew = xk + alpha * pk;
    fnew = f(xnew);
    
    c1_gradfk_pk = c1 * gradfk' * pk;
    bt = 0;

    % Backtracking strategy: 
    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)

        alpha = rho * alpha;

        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        bt = bt + 1;
    end

    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        break
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    if nargin == 9
        gradfk = gradf(xk);
    else
        gradfk = gradf(xk, hstep, hstep_i);
    end

    gradfk_norm = norm(gradfk);
    
    k = k + 1;
    
    xseq(:, k) = xk;

    btseq(k) = bt;
end

% Cut xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);

% Add x0 at the beginning of xseq
xseq = [x0, xseq];

end