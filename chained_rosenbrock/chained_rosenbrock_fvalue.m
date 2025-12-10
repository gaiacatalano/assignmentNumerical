function f = chained_rosenbrock_fvalue(x)
    
    x_im1 = x(1:end-1);   % x(i-1)
    x_i   = x(2:end);     % x(i)

    f = sum( 100*(x_im1.^2 - x_i).^2 + (x_im1 - 1).^2 );


    % n = length(x);
    % 
    % f = 0;
    % 
    % for i = 2:n
    % 
    %     % xim1 = x(i-1);
    %     % xi = x(i);
    % 
    %     f = f + 100*(x(i-1)^2 - x(i))^2 + (x(i-1) - 1)^2;
    % 
    % end

end