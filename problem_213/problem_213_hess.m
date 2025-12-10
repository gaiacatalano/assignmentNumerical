% Computes the Hessian of the Problem 213

function hess = problem_213_hess(x)
    
    n = length(x);
    h = 1/(n+1);

    x = [0; x; 1];
    xkm1 = x(1:n);
    xk   = x(2:n+1);
    xkp1 = x(3:n+2);

    fk = 2*xk + h^2*(xk + sin(xk)) - xkm1 - xkp1;
    
    % derivate df_k/dx_k e d^2 f_k/dx_k^2
    df_dxk   = 2 + h^2*(1 + cos(xk));   % A_k
    d2f_dxk2 = -h^2 * sin(xk);          % B_k

    % diagonale principale d0
    % contributo da k=j: df^2 + f * d2f
    d0 = fk .* d2f_dxk2 + df_dxk.^2;

    % contributo da f_{j-1} e f_{j+1} (df/dx_j = -1)
    if n > 1
        d0(1:n-1) = d0(1:n-1) + 1;  % da f_{j+1}
        d0(2:n)   = d0(2:n)   + 1;  % da f_{j-1}
    end
    
    % prima sotto-diagonale dm1
    % H_{j,j+1} = -(df_dxk(j) + df_dxk(j+1))
    dm1 = zeros(n,1);
    dm1(1:n-1) = -(df_dxk(1:n-1) + df_dxk(2:n));

    % seconda sotto-diagonale dm2
    % H_{j,j+2} = 1 (da f_{j+1})
    dm2 = zeros(n,1);
    dm2(1:n-2) = 1;
    
    hess = spdiags([dm2 dm1 d0 [0; dm1(1:end-1)] [0; 0; dm2(1:end-2)]], ...
                   [-2  -1   0   1                2], n, n);

    % simmetrizzazione numerica
    hess = 0.5 * (hess + hess.');

end

%     d0 = zeros(n,1);
%     dp1 = zeros(n,1);
%     dp2 = zeros(n,1);
%     dm1 = zeros(n,1);
%     dm2 = zeros(n,1);
% 
%     for i = 1:n
% 
%         fk = 2*x(i+1) + (h^2)*(x(i+1) + sin(x(i+1)))-x(i)-x(i+2);
%         d0(i) = (2 + (1 + cos(x(i+1))*h^2))^2 - sin(x(i+1))*h^2*fk;
% 
%         if i ~= 1
%             d0(i) = d0(i) + 1;
%         end
% 
%         if i ~= n
%             d0(i) = d0(i) + 1;
%             dm1(i) = - 4 - (2 + cos(x(i+1)) + cos(x(i+2)))*h^2;
%         end
% 
%         if i < n-1
%             dm2(i) = 1;
%         end
% 
%     end
% 
%     hess = spdiags([dm2 dm1 d0 [0; dp1(1:end-1)] [0; 0; dp2(1:end-2)]], [-2 -1 0 1 2], n, n);
%     hess = 0.5 * (hess + hess');
% 
% end