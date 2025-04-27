function[F] = discrete_boundary_value_problem(n, x_bar, h)
    
    x_first = 0;
    x_last = 0;
    x = [x_first; x_bar; x_last];

    f = zeros(n,1);
    for i=1:n
        f(i)= (2*x(i+1)-x(i)-x(i+2) + (h^2 * (x(i+1) + i*h + 1)^3)/2)^2;
    end

    F = @(x) sum(f.^2);

end