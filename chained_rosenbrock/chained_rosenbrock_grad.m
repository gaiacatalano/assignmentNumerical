function grad = chained_rosenbrock_grad(x)
    
    n = length(x);
    grad = zeros(n,1);

    x_im1 = x(1:end-1);   % x(1),...,x(n-1)
    x_i   = x(2:end);     % x(2),...,x(n)

    % contributo quando la variabile appare come x_{i-1} nel termine i
    % per i = 2,...,n  → indici 1,...,n-1 in x
    left = 400 * x_im1 .* (x_im1.^2 - x_i) + 2 * (x_im1 - 1);

    % contributo quando la variabile appare come x_i nel termine i
    % per i = 2,...,n → indici 2,...,n in x
    right = -200 * (x_im1.^2 - x_i);

    % sommo i contributi nei rispettivi indici
    grad(1:end-1) = grad(1:end-1) + left;    % include il caso i = 1
    grad(2:end)   = grad(2:end)   + right;   % include il caso i = n
    
end

%     n = length(x);
% 
%     grad = zeros(n,1);
% 
%     for i = 2:n-1
%         % xim1 = x(i-1);
%         % xi = x(i);
%         % 
%         % df_dxim1 = 400*(xim1^2 - xi)*xim1 + 2*(xim1 - 1);
%         % df_dxi   = -200*(xim1^2 - xi);
% 
%         % df_dxim1 = 400*(xi^2 - x(i+1))*xi + 2*(xi - 1);
%         % df_dxi   = -200*(xim1^2 - xi);
%         % 
%         % grad(i-1) = grad(i-1) + df_dxim1;
%         % grad(i) = grad(i) + df_dxi;
%         grad(i) = 400*x(i)*(x(i)^2 - x(i+1)) + 2*(x(i) - 1) - 200*(x(i-1)^2 - x(i));
% 
% 
%     end
% 
%     grad(1) =  400*x(1)*(x(1)^2 - x(2)) + 2*(x(1) - 1);
%     grad(n) =  -200*(x(n-1)^2 - x(n));
% 
% end