function [s, z, Bx] = SimplexMatricial(A, c, b, Bx)
    [m,n] = size(A);
    tol = 1e-8; % Tolerancia mayor a eps para estabilidad
    % Primer factorización para entrar al while
    B = A(:, Bx);
    [L, U, p] = lups(B);
    % Se calcula lam para obtener r
    cB = c(Bx);
    w = forsub(U', cB);
    v = backsub(L', w);         % Se calculan los multiplicadores
    lambda = zeros(m, 1); 
    lambda(p) = v;              % Se aplica permutación inversa
    
    r = c - A'*lambda;
    noBas = setdiff(1:n, Bx); % Índices de las variables no básicas
    [R, idx] = max(r(noBas));
    vin = noBas(idx);
    
    aj = A(:, vin);
    d = backsub(U, forsub(L, aj(p)));

    condition = true;
    while condition
        % Se calcula xB para el test de razón
        xB = backsub(U, forsub(L, b(p)));
        % Test de razón mínima
        razones = xB ./ d;
        razones(d <= tol) = inf; 
        [~, vout] = min(razones);      
        % Actualizar la base
        Bx(vout) = vin;
        
        % Se actualiza L
        B = A(:, Bx);
        [L, U, p] = lups(B);
        
        % Recalcular lam y r
        cB = c(Bx);
        w = forsub(U', cB);
        v = backsub(L', w);
        lambda = zeros(m, 1); lambda(p) = v;
        r = c - A'*lambda;
        noBas = setdiff(1:n, Bx);
        [R, idx] = max(r(noBas));
        vin = noBas(idx);
        
        % Actualizar d para la nueva variable de entrada y condición
        aj = A(:, vin);
        d = backsub(U, forsub(L, aj(p)));
        condition = (R > tol) && any(d > tol);
    end

    % Manejo de acotamiento y múltiples soluciones
    if (R > tol) && all(d <= tol)
        error('El problema NO ESTÁ ACOTADO');
    end
    if any(abs(r(noBas)) < tol)
        disp('Se han detectado múltiples soluciones óptimas');
    end

    % Cálculo final con la base óptima encontrada
    xB = backsub(U, forsub(L, b(p)));
    s = zeros(n, 1);
    s(Bx) = xB;
    z = c' * s;
    if z < 0
        z=-z;
    end
end