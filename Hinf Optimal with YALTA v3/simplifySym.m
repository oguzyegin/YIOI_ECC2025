function [L,Lnum,Lden]=simplifySym(L)
syms s
L = simplify(L);
[Lnum,Lden]=numden(L);
Lnum=sym2poly(Lnum);
Lden=sym2poly(Lden);
Lnum=round(Lnum/Lden(1)*1e+4)*1e-4;
Lden=round(Lden/Lden(1)*1e+4)*1e-4;
Lnum=poly2sym(Lnum,s);
Lden=poly2sym(Lden,s);
L = Lnum/Lden;