function TF = sym2zpk(X)
[nX,dX] = numden(X);
nX = sym2poly(nX);
dX = sym2poly(dX);
TF = zpk(roots(nX),roots(dX),nX(1)/dX(1),'DisplayFormat','Roots');