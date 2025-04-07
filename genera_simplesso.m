function [simplex_points, y] = genera_simplesso(x0)
    % x0: punto iniziale 2x1
    % simplex_points: matrice 2x3, ogni colonna è un vertice del simplesso
    % y: punto generato all'interno del simplesso

    % Parametro di perturbazione
    delta = 0.1;

    % Genera gli altri due punti del simplesso
    x1 = x0 + [delta; 0];
    x2 = x0 + [0; delta];

    % Costruisci il simplesso: ogni colonna è un punto in R^2
    simplex_points = [x0, x1, x2];

    % Genera lambda casuali tali che somma = 1 e lambda_i >= 0
    lambda = rand(3,1);
    lambda = lambda / sum(lambda);  % Normalizza

    % Calcola il punto interno y
    y = lambda(1) * simplex_points(:,1) + ...
        lambda(2) * simplex_points(:,2) + ...
        lambda(3) * simplex_points(:,3);
end
