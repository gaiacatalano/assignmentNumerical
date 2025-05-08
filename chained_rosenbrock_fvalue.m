function f = chained_rosenbrock_fvalue(x)

    n = length(x);

    f = 0;

    for i = 2:n

        xim1 = x(i-1);
        xi = x(i);

        f = f + 100*(xim1^2 - xi)^2 + (xim1 - 1)^2;

    end

end