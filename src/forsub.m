function y = forSub(L, bp)
    [n,~]=size(L);
    y = zeros(n,1);
    y(1) = bp(1)/L(1,1);
    for i = 2 : n
        y(i) = (bp(i) - L(i, 1:i-1) * y(1:i-1)) / L(i,i);
    end
end