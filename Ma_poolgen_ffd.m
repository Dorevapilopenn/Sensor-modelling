function r=Ma_poolgen(P)
[a,b,c,d,e,f,g,h,i,j] = ndgrid(P,P,P,P,P,P,P,P,P,P);
f=[a(:), b(:), c(:), d(:), e(:), f(:), g(:), h(:), i(:), j(:)]; 
r = f(f(:,3) >= f(:,2), :);
end