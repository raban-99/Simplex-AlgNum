function [s,z,Bx] = SimplexMatricial(A,c,b,Bx)
    [~,n] = size(A);
    r = inf*ones(n,1);     % r = costos reducidos
    d = 1;                 % d = direccion simplex
    theta = inf;           % theta = razon minima
    condition = true;

    while condition
        B = A(:,Bx);
        % Factorizacion LU de la base
        [L,U,p,~] = lups(B);

        % Resolver B*xB = b
        bp = b(p);
        y = forsub(L,bp);
        xB = backsub(U,y);

        % Vector solucion completo
        s = zeros(n,1);          % s = solución del problema
        s(Bx) = xB;

        % Costos basicos
        cB = c(Bx);              % cB = costos de las variables basicas

        % Resolver B' * lam = cB
        w = backsub(U',cB);
        lam = forsub(L',w);      % lam = multiplicadores simplex
        % Costos reducidos
        r = c - A'*lam;          % r = costos reducidos

        % Variable que entra a la base
        [~,input] = max(r);      % input = índice de variable entrante

        % Dirección simplex  B*d = aj
        aj = A(:,input);         % aj = columna a analizar
        y = forsub(L,aj(p));
        d = backsub(U,y);        % d = dirección simplex

        % Test de razón mínima
        raz = xB ./ d;           % raz = razones para elegir variable que sale
        raz(d <= eps) = inf;
        [theta,output] = min(raz);   % output = índice de variable de salida

        % Actualizar base y condición
        Bx(output) = input;
        condition = max(r) > eps && theta < inf;
    end

    % Detectar problema no acotado
    if theta == inf
        error('El problema no está acotado')
    end
    % Detectar múltiples soluciones óptimas
    noBas = setdiff(1:n,Bx);     % variables no básicas
    if any(abs(r(noBas)) <= eps)
        disp('Existen múltiples soluciones óptimas')
    end

    % Valor de la funcion objetivo
    z = c'*s;
end