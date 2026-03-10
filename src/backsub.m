function x = backsub(A,b)
n = length(b);
x = zeros(n,1);
for i = n:-1:1
    suma = 0;
    for j = i+1:n
        suma = suma + A(i,j)*x(j);
    end
    x(i) = (b(i) - suma)/A(i,i);    
end
end