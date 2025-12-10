function hess = chained_rosenbrock_hess(x)

    n = length(x);
    hess = spalloc(n,n,3*n-2);

    x_im1 = x(1:end-1);   % x(1),...,x(n-1)
    x_i   = x(2:end);     % x(2),...,x(n)

    d0 = zeros(n,1);
    %dm1 = zeros(n,1);
    dm1 = -400*x_im1;

    % contributi g_aa su x_{i-1} (i = 2..n) → indici 1..n-1
    d0(1:end-1) = d0(1:end-1) + 400*(3*x_im1.^2 - x_i) + 2;

    % contributi g_bb su x_i (i = 2..n) → indici 2..n
    d0(2:end)   = d0(2:end)   + 200;

    hess = spdiags([[0; dm1], d0, [dm1;0]], [-1 0 1], n, n);
    hess = 0.5*(hess + hess');

    % for i = 2:n
    % 
    %     % if i ~= 1
    %     %     d0(i) = d0(i) + 200;
    %     %     dm1(i) = -400*x(i-1);
    %     % end
    % 
    %     if i ~= n
    %         %d0(i) = d0(i) + 400*(3*x(i)^2 - x(i+1)) + 2;
    %         d0(i) = 400*(3*x(i)^2 - x(i+1)) + 202;
    %     end
    %     dm1(i-1) = -400*x(i-1);
    % 
    %     % xim1 = x(i-1);
    %     % xi = x(i);
    %     % 
    %     % d2f_dxim12 = 400*(3*xim1^2 - xi) + 2;
    %     % d2f_dxim1xi = -400*xim1;
    %     % d2f_dxi2 = 200;
    % 
    %     % d0(i-1) = d0(i-1) + d2f_dxim12;
    %     % d0(i)   = d0(i)   + d2f_dxi2;
    %     % dp1(i) = d2f_dxim1xi;
    %     % dm1(i-1) = d2f_dxim1xi;
    % 
    % 
    % end
    % 
    % d0(1) = 400*(3*x(1)^2 - x(2)) + 2;
    % d0(n) = 200;
    % 
    % hess = spdiags([dm1, d0, [0; dm1(1:end-1)]], [-1 0 1], n, n);
    % hess = 0.5*(hess + hess');
end