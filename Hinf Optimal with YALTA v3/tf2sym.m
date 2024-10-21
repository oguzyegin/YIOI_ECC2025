function TF = tf2sym(X)
syms s;
n = poly2sym(X.num{1},s);
d = poly2sym(X.den{1},s);
TF = n/d;