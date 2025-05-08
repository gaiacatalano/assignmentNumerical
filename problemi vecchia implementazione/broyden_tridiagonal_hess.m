function hess = broyden_tridiagonal_hess(x)

    n = length(x);

    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];

    d0 = zeros(n,1);
    dp1 = zeros(n,1);
    dp2 = zeros(n,1);

    for k=1:n

        xkm1 = x(k);
        xk = x(k+1);
        xkp1 = x(k+2);
        
        df_dxkm1 = -1;
        df_dxk   = 3 - 4*xk;
        df_dxkp1 = -2;

        fk = (3-2*xk)*xk - xkm1 -2*xkp1 + 1;

        % Diagonale principale
        if k > 1
            d0(k) = d0(k) + df_dxkp1^2;
        end
        d0(k) = d0(k) + df_dxk^2 + fk * (-4);
        if k < n
            d0(k) = d0(k) + df_dxkm1^2;
        end

        % Bande fuori diagonale
        if k < n
            dp1(k) = dp1(k) + df_dxk * df_dxkp1 + df_dxkm1 * df_dxk;
        end
        if k <= n-2
            dp2(k) = dp2(k) + df_dxkp1^2;
        end
    end

    hess = spdiags([dp2 dp1 d0 dp1 dp2], [-2 -1 0 1 2], n, n);
    hess = 0.5 * (hess + hess');

end
