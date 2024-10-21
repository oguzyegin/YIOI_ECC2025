function K_opt = findKopt(beta,gamma_opt,Ro,Mn,M1, G_hat)
% beta = 8.903931715820827j;
% G_hat = hat_G_gamma;
[~,den_Ro] = numden(Ro);
poles_Ro = roots(sym2poly(den_Ro));
uns_poles_Ro = poles_Ro(real(poles_Ro)>0);
uns_poles_Ro = [uns_poles_Ro; beta];
nu = length(uns_poles_Ro);
results = zeros(1,nu);
rts = zeros(1,nu);
A = [];
for k=1:nu
    s = uns_poles_Ro(k);
    results(k) = -eval((Mn*M1*G_hat)*gamma_opt);
    v1 = (-1).^(1:nu)-results(k);
    v2 = s.^((nu-1):-1:0);
    B(k,1) = (results(k)+1)*s^(nu);
    A = [A; v1.*v2];
end
K_opt = A\B