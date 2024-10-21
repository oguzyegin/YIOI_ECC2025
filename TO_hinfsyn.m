function [Copt, g, Ca, Hd, Hda] = TO_hinfsyn(nw1,dw1,nw2,dw2,a,Mn,Md,No,Hda_ord)
close all
thr = 1e-8;
g = gammaOpt(nw1,dw1,nw2,dw2,a,1e3,...
    0.1,10,thr,Mn);
%%
if nargin<9
    order_inputs.Hd_order = 7;
else
    order_inputs.Hd_order = Hda_ord;
end
order_inputs.Hd_rel_order = 1;
order_inputs.DO_APP_EVEN_FINITE=0;
lw1 = -3; lw2 = 4; Lpoints = 1e3;
[warning, E_gamma, G_gamma, F_gamma, L_gamma, Ro, Hn, Hda, Hd, hat_G_gamma, per, Ca, app_No, Copt] = optController(g,nw1,dw1,nw2,dw2,a,thr,lw1,lw2,Lpoints,Mn,Md,No,order_inputs);