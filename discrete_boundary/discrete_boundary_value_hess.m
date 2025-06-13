function hess = discrete_boundary_value_hess(x)
    
    n = length(x);
    h = 1/(n+1);

    x_first = 0;
    x_last = 0;
    x = [x_first; x; x_last];    

    d0 = zeros(n,1);
    %dp1 = zeros(n,1);
    %dp2 = zeros(n,1);
    dm1 = zeros(n,1);
    dm2 = zeros(n,1);

    for k=1:n

        % xkm1 = x(k);
        % xk = x(k+1);
        % xkp1 = x(k+2);

        fk = (2*x(k+1)-x(k)-x(k+2) + (h^2 * (x(k+1) + k*h + 1)^3)/2)^2;
        d0(k) = 2*(2 + 3/2 * h^2 * (x(k+1) + k*h + 1)^2)^2 + 2*fk*(3*h^2 * (x(k+1) + k*h + 1));

        % df_dxkm1 = -1;
        % df_dxk   = 2 + 3/2 * h^2 * (xk + k*h + 1)^2;
        % df_dxkp1 = -1;

        % d2f_dxdk2 = 3*h^2 * (xk + k*h + 1);

        if k ~= 1
            d0(k) = d0(k) + 1;            
        end
        
        if k ~= n
            d0(k) = d0(k) + 1;
            dm1(k) = -(2 + 3/2 * h^2 * (x(k+1) + k*h + 1)^2) - (2 + 3/2 * h^2 * (x(k+2) + k*h + 1)^2);
        end

        if k < n-1
            dm2(k) = 2;
        end

        
        % if k < n
        %     dm1(k) = dm1(k) + 2*(df_dxk*df_dxkp1 + df_dxkm1*df_dxk);
        %     dp1(k+1) = dp1(k+1) + 2*(df_dxk*df_dxkp1 + df_dxkm1*df_dxk);
        % end
        % if k < n-1
        %     dm2(k) = dm2(k) + 2*df_dxkp1^2;
        %     dp2(k+2) = dp2(k+2) + 2*df_dxkp1^2;
        % end

    end

    %hess = spdiags([dm2 dm1 d0 dp1 dp2], [-2 -1 0 1 2], n, n);
    hess = spdiags([dm2 dm1 d0 [0; dm1(1:end-1)] [0; 0; dm2(1:end-2)]], [-2 -1 0 1 2], n, n);
    hess = 0.5 * (hess + hess');
    %condhess = condest(hess)

end