function [Hn,Hd]=findHnHd(gamma_opt,K_opt,Ro,Mn,M1,G_hat,dInf)
syms s;
[num,den]=numden(K_opt);
num = sym2poly(num);
den = sym2poly(den);
[~,polesK,~]=residue(num,den);
[num2,den2]=numden(K_opt*Ro/gamma_opt/dInf);
den2 = sym2poly(den2);
num2 = sym2poly(num2)/den2(1);
den2 = den2/den2(1);
[res,poles,~]=residue(num2,den2);
polesHn = polesK(real(polesK)>0);

[GS,GNS]=stabsep(tf(num2,den2));

Hn = 0*s;
if ~isempty(polesHn)
    for k=1:length(polesHn)
        a = res(abs(polesHn(k)-poles)<1e-3);
        newH = a(1)/(s-polesHn(k));
        Hn = Hn + newH;
    end
end
[num_Hn,den_Hn] = numden(Hn);
dHn = real(sym2poly(den_Hn));
nHn = real(sym2poly(num_Hn))/dHn(1);
dHn = dHn/dHn(1);
Hn = poly2sym(nHn,s)/poly2sym(dHn,s);
decompositionEqn = (Ro/dInf)*(K_opt/gamma_opt+Mn*M1*G_hat)-1;
Hd = decompositionEqn-Hn;