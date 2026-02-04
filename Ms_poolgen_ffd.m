function r=Ms_poolgen(P)
[A,B,C,D,E] = ndgrid(P,P,P,P,P);
f=[A(:), B(:), C(:), D(:), E(:)];
r = f(f(:,3) >= f(:,2), :);
end