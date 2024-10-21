clear all
close all
clc
addpath 'Hinf Optimal with YALTA v3';
addpath 'IRKA';
%% Optimal Controller & Perf Level for 
aa = 1; khat = 2; hhat = 0.75;

W1 = tf(1,[1 0]);
W2 = tf([khat 0],1);

syms s;
Mn = exp(-hhat*s);
Md = (s-aa)/(s+aa);
No = 1/(s+aa);

[Copt, g, ~, Hd, ~] = TO_hinfsyn(W1.num{1},W1.den{1},...
    W2.num{1},W2.den{1},aa,Mn,Md,No);

gopt = g;
kghat=khat*sqrt((2/khat)-gopt^(-2));
zg = g*exp(-hhat)/(khat+kghat+1);
phopt = (zg-1)/(zg+1);

% Hopt
syms s;
Hopt = (s*(phopt*s-gopt^(-2))+(s^2+gopt^(-2))+...
  (1/khat)*( (s+phopt)*(kghat*s+1)+gopt*s*(phopt-s)*exp(-hhat*s)))/...
       ( (s^2+gopt^(-2))*(s-1) );
sx = tf('s');
HoptTF= (sx*(phopt*sx-gopt^(-2))+(sx^2+gopt^(-2))+...
      (1/khat)*( (sx+phopt)*(kghat*sx+1)+gopt*sx*(phopt-sx)*exp(-hhat*sx)))/...
           ( (sx^2+gopt^(-2))*(sx-1) );

% PI block of opt. controller
Cpi = gopt*tf([1 -phopt],[khat 0]);

w = logspace(-3,4,1e3); s = 1j*w;
Hopt_w = eval(Hopt);

desired_order = 6;
% Pade approximation
[np,dp] = pade(hhat,desired_order);
Hpade = minreal((sx*(phopt*sx-gopt^(-2))+(sx^2+gopt^(-2))+...
      (1/khat)*( (sx+phopt)*(kghat*sx+1)+gopt*sx*(phopt-sx)*tf(np,dp)))/...
           ( (sx^2+gopt^(-2))*(sx-1) ),1e-3);
% tfest 40 order
opt = tfestOptions('EnforceStability',1);
Hopt_tfest = tfest(frd(Hopt_w,w),40,opt);

% SSEST
Hopt_ssest = ssest(frd(Hopt_w,w),desired_order);

%% IRKA
F = ss(Hopt_tfest);
close all;
nred = desired_order;
F_h = hankelmr(Hopt_tfest,nred);
[SV,w] = sigma(F-F_h);
[pks,locs] = findpeaks(SV);
[~,idx] = sort(pks,'descend');
b = randn(size(F.B,2),nred);
c = randn(nred,size(F.C,1));

[V,W,A,B,C,E,s_interp]=IRKA(F.A,... % system Mat. FOM ex.
    F.B,...                                  % input  Mat. FOM ex.
    F.C,...                                  % output Mat. FOM ex.
    eye(size(F.A,1)),...                   % descriptor Matrix
    sort(w(locs(idx(1:nred))),'ascend'),...       % initial shift value
    b,...                                         % tangential direction for B
    c,...                                         % tangential direction for C
    25,...                                        % maximum iteration number
    1.0e-3,...                                    % convergence tolerance
    1);                                           % debug flag
clear Vr Wr;
N = size(F.A,1);
Vr = V;
Wr = W;
Er = E;
Ar = A;
br = B;
cr = C';

Hr0 = ss(Ar,br,cr',[]); Hr0.E = Er;
%%

% BALRED
opts = balredOptions('StateProjection','Truncate');
H_balred_tfest = balred(prescale(ss(Hopt_tfest)),desired_order,opts);

% MM results
syms s;
H_m6 = (0.5732*s^5 + 137.8*s^4 - 384.1*s^3 + 5.263e04*s^2 - 1.039e05*s + 3.044e06)/...
    (s^6 + 33.7*s^5 + 750.2*s^4 + 1.042e04*s^3 + 9.543e04*s^2 + 5.248e05*s + 1.335e06);
H_m6 = sym2tf(H_m6);
%%
w = logspace(-2,3.3,1e4);

H_hna = tf(hankelmr(Hopt_tfest,desired_order));

close all;
C_ssest = Cpi*feedback(1,Hopt_ssest);
[psi_ssest,g_ssest] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_ssest,w);

C_balred = Cpi*feedback(1,H_balred_tfest);
[psi_balred,g_balred] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_balred,w);

C_hna = Cpi*feedback(1,H_hna);
[psi_hna,g_hna] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_hna,w);

% SYSTUNE
sys = ss(Hopt_tfest);
ny = size(sys,1); nu = size(sys,2);
nr= desired_order;
opt1 = systuneOptions('SoftTol',1e-6);
K1 = tunableSS('Controller',nr,ny,nu); rsys_bt = balred(sys,nr);
K1.a.Value = rsys_bt.a; K1.b.Value = rsys_bt.b; K1.c.Value = rsys_bt.c; 
K1.d.Value = 0; K1.D.Maximum = 0;
cl1 = sys-K1;
cl1.InputName = 'w';
cl1.OutputName = 'z';
HinfReq = TuningGoal.WeightedGain('w','z',1,tf([1 1],[1/50 1 1e-1]));
opt1.Hidden.Problem = 'Hinf';
[cl2,f2] = systune(cl1,HinfReq,[],opt1);

rsys_systune = ss(cl2.Blocks.Controller);

C_systune = Cpi*feedback(1,rsys_systune);
[psi_systune,g_systune] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_systune,w);

C_irka = Cpi*feedback(1,Hr0);
[psi_irka,g_irka] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_irka,w);

C_pade = Cpi*feedback(1,Hpade);
[psi_pade,g_pade] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_pade,w);

C_m6 = Cpi*feedback(1,H_m6);
[psi_mm6,g_mm6] = perf_level(W1,W2,sym2tf(No/Md),Mn,C_m6,w);

% figure(5); clf; semilogx(w,psi_ssest,'c','LineWidth',2);
% hold on; semilogx(w,psi_balred, 'b','LineWidth',2);
% semilogx(w,psi_hna,'k','LineWidth',2);
% semilogx(w,psi_systune,'g','LineWidth',2);
% semilogx(w,psi_irka,'k-.','LineWidth',2);
% semilogx(w,psi_pade,'b-.','LineWidth',2);
% semilogx(w,psi_mm6,'r','LineWidth',2);
% legend({'ssest','tfest-balred','hna','systune','IRKA','pade','moment matching'},'Location','southeast')
% xlabel('Frequency (rad/s)'); ylabel('$\Psi_a(j\omega)$','Interpreter','latex');
% title('Performance Levels for $\nu=6,~h=0.75$','Interpreter','latex')
% set(gca,'FontSize',16,'FontName','Times New Roman'); grid on;

%% Plant approximation
[np,dp] = pade(hhat,desired_order);
Mn_a = tf2sym(tf(np,dp));
Delta_m = ((Mn-Mn_a))/Mn;
s = 1j*w; m_Dm = 20*log10(abs(eval(Delta_m))); 
figure(6); clf; semilogx(w,m_Dm); hold on;
W_m1 = tf(0.0057*[1 0 0],[1/20^2 1.2/20 1]); [m_m1,~] = bode(W_m1,w); m_m1 = 20*log10(squeeze(m_m1));
semilogx(w,m_m1);
W1 = tf(1,[1 0]);
W2 = tf([khat 0],1);
Wma = W2+W_m1;

[Ca, ga] = TO_hinfsyn(W1.num{1},W1.den{1},Wma.num{1},Wma.den{1},1,Mn_a,Md,No,7);
Ca = sym2tf(Ca);
%% Plot all
Ca = minreal(Ca,1e-2);
[psi_pa,g_pa] = perf_level(W1,W2,sym2tf(No/Md),Mn,Ca,w);
close all;

% figure(5); semilogx(w,psi_ssest,'c','LineWidth',2); hold on; 
% figure(5); semilogx(w,psi_balred, 'b','LineWidth',1); hold on; 
% figure(5); hold on; semilogx(w,psi_hna,'k','LineWidth',2);
% figure(5); hold on; semilogx(w,psi_pade,'k-.','LineWidth',2);
% figure(5); hold on; semilogx(w,psi_pa,'b-.','LineWidth',2);
% legend({'ssest','tfest-balred','hna','systune','IRKA','pade','plant app.','moment matching'},'Location','northwest')
figure(5); semilogx(w,psi_systune,'k','LineWidth',1);
hold on; semilogx(w,psi_irka,'g','LineWidth',1);
semilogx(w,psi_mm6,'r','LineWidth',1);
semilogx(w([1 length(w)]),[gopt gopt],'b');
legend({'systune','IRKA','MM','\gamma_{opt}'},'Location','northwest')
xlabel('Frequency (rad/s)'); ylabel('$\Psi_a(j\omega)$','Interpreter','latex');
title('Performance Levels for $\nu=6,~h=0.75$','Interpreter','latex')
set(gca,'FontSize',16,'FontName','Times New Roman'); grid on;

g6_h75 = [g_ssest, g_balred, g_hna, g_systune, g_irka, g_pade, g_mm6, g_pa];
xlim([1e-1,1e3]); ylim([7,11.5]);
%%
Cpi_pa = Ca.num{1}(1)*tf([1, 0.08144],[1 0]);
H_temp = minreal(Cpi_pa/Ca,1e-3);
H_pa = tf(H_temp.num{1}-H_temp.den{1},H_temp.den{1});
zpk(H_pa)

%%
H0 = feedback(1,HoptTF);
Copt_TF = H0*Cpi;
[m1,p1] = bode((Copt_TF-C_systune)/Copt_TF,w); m1=squeeze(m1); p1=squeeze(p1);
[m2,p2] = bode((Copt_TF-C_m6)/Copt_TF,w); m2=squeeze(m2); p2=squeeze(p2);
[m3,p3] = bode((Copt_TF-C_irka)/Copt_TF,w); m3=squeeze(m3); p3=squeeze(p3);
figure(6); clf; loglog(w,m1,'k','LineWidth',2);
hold on; loglog(w,m3,'g','LineWidth',2); loglog(w,m2,'r','LineWidth',2);
legend({'systune','IRKA','MM'},'Location','northwest')
xlabel('Frequency (rad/s)'); ylabel('Relative Error','Interpreter','latex');
title('$|| 1 - C_a/C_{opt}||_\infty$ for $\nu=6,~h=0.75$','Interpreter','latex')
set(gca,'FontSize',16,'FontName','Times New Roman'); grid on;
xlim([1e-1,1e3]); ylim([5e-5,1]);

