function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    modified_newton_bcktrck_preconditioning(x0, f, gradf, Hessf, ...
    kmax, tolgrad, c1, rho, btmax)

% Function that performs the Newton optimization method, using
% backtracking strategy for the step-length selection.

% Function handle for the armijo condition
farmijo = @(fk, alpha, c1_gradfk_pk) ...
    fk + alpha * c1_gradfk_pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
gradfk = gradf(xk);
k = 0;
gradfk_norm = norm(gradfk);
Hessfk = Hessf(xk);

delta = sqrt(eps);


while k < kmax && gradfk_norm >= tolgrad

    if issparse(Hessfk)
        if any(isnan(Hessfk(:))) || any(isinf(Hessfk(:)))
            error('Hessiana contiene NaN o Inf alla iterazione %d', k);
        end
        lambda_min = eigs(Hessfk, 1, 'SA');
    else
        lambda_min = min(eig(Hessfk));
    end
    tau_k = max(0, delta-lambda_min);
    E_k = tau_k*speye(size(Hessfk,1));
    B_k = Hessfk + E_k;
    
    % % solved with pcg with preconditioning, con B_k=LL' incomplete Choleski
    % L = ichol(B_k);
    % [pk, ~, ~, iterk, ~] = pcg(B_k, -gradfk, [], [], L, L');

    % solved with pcg with diagonal preconditioning
    M_diag = diag(B_k);
    M_diag(M_diag <= 0) = 1;  % evita divisioni per zero o valori negativi
    n = length(xk);
    M = spdiags(M_diag, 0, n, n);
    [pk, ~] = pcg(B_k, -gradfk, 1e-6, kmax, M);

    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    c1_gradfk_pk = c1 * gradfk' * pk;

    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
    end
    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        break
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
% "Add" x0 at the beginning of xseq (otherwise the first el. is x1)
xseq = [x0, xseq];

end