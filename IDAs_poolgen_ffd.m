function r=IDAs_poolgen(P)
[A,B,C] = ndgrid(P,P,P);
f=[A(:), B(:), C(:)];
r = f(f(:,3) >= f(:,2), :);
end