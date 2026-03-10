function [L,U,p,s] = lups(A)
    [n,~] = size(A);
    % Asegurarse de que es una matriz cuadrada error
    L = eye(n);
    p = 1:n;
    s = 1;
    for j = 1 : n-1
        [~,r] = max(abs(A(j:n,j))); %Verificar que no haya ceros sino inf sol.
        mayor = j + (r - 1); % Ajustar el índice para el pivoteo
        if abs(A(mayor,j)) < eps %Verificar |pivote|<eps error matriz singular
            error('Matriz singular');
        end
        if j ~= mayor
            A([j, mayor], :) = A([mayor, j], :);
            L([j, mayor], 1:j-1) = L([mayor, j], 1:j-1); % Se actualizan los índices de los pivotes
            p([j,mayor]) = p([mayor, j]); % Se guardan los multiplicadores en L
            s=-s;
        end
        for i = j+1 : n
            L(i,j) = A(i,j)/A(j,j);
            A(i,j:n)=A(i,j:n) - L(i,j)*A(j,j:n);
        end
    end    
    U=A; % Se copia una sola vez en este punto
end