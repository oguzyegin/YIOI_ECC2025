function TF = sym2tf(X)
[nX,dX] = numden(X);
dX = sym2poly(dX);
nX = sym2poly(nX)/dX(1);
dX = dX/dX(1);
TF = tf(nX,dX);