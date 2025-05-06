function [f, gradf, Hess] = chained_rosenbrock(n, x_bar)

    f = zeros(n,1);
    grad = zeros(n, 1);

    d0 = zeros(n, 1);   
    dp1 = zeros(n-1, 1); 
    dm1 = zeros(n-1, 1); 

    for i = 2:n
        xim1 = x_bar(i-1);
        xi = x_bar(i);

        f(i) = 100*(xim1^2 - xi)^2 + (xim1 - 1)^2;

        % Gradient
        df_dxim1 = 400 * (xim1^2 - xi) * xim1 + 2 * (xim1 - 1);
        df_dxi   = -200 * (xim1^2 - xi);

        grad(i-1) = grad(i-1) + df_dxim1;
        grad(i)   = grad(i)   + df_dxi;

        % Hessiana

        d2f_dxim12 = 400 * (3 * xim1^2 - xi) + 2;
        d2f_dxim1xi   = -400 * xim1;
        d2f_dxi2     = 200;

        d0(i-1) = d0(i-1) + d2f_dxim12;
        d0(i)   = d0(i)   + d2f_dxi2;
        dp1(i-1) = d2f_dxim1xi;
        dm1(i-1) = d2f_dxim1xi;

    end
    % usando H = spdiags([dm1 d0 dp1], [-1 0 1], n, n); mi dava errore
    % perche dm1, d0 e dp1 hanno lunghezza n-1, n e n-1, quindi non posso
    % concatenare perche hanno lunghezze diverse
    % dice di riempire dm1 e dp1 con zero per arrivare ad n

    % Allunga dm1 e dp1 a lunghezza n
    dm1_full = [dm1; 0];
    dp1_full = [0; dp1];
    
    % Costruzione Hessiana sparsa
    H = spdiags([dm1_full d0 dp1_full], [-1 0 1], n, n);
    
    F = @(x) sum(f);
    gradf = @(x) grad;
    Hess = @(x) H;

end
