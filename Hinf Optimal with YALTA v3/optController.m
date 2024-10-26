% % Computes optimal and suboptimal controller by using new structure
function [warning, E_gamma, G_gamma, F_gamma, L_gamma, Ro, Hn, app_Hd, Hd, hat_G_gamma, per, Ca, app_No, Copt]=optController(gamma_opt,nw1,dw1,nw2,dw2,alpha,eps,lw1,lw2,Lpoints,Mn,Md,No,order_inputs)

syms s;
Ro = nan;
Hn = nan;
EPSSS=eps;
Plant = Mn*No/Md;
% Calculation of W1(s)
W1_n = poly2sym(nw1,s);
W1_d = poly2sym(dw1,s);
W1_s = W1_n/W1_d;
% Calculation of W2(s)
W2_n = poly2sym(nw2,s);
W2_d = poly2sym(dw2,s);
W2_s = W2_n/W2_d;

nwm1=zeros(1,length(nw1));
dwm1=zeros(1,length(dw1));
nwm2=zeros(1,length(nw2));
dwm2=zeros(1,length(dw2));
% Calculation of W1(-s)
for k = 1:length(nw1)
    nwm1(k) = nw1(k)*(-1)^(length(nw1)-k);
end
for k = 1:length(dw1)
    dwm1(k) = dw1(k)*(-1)^(length(dw1)-k);
end
mW1_n = poly2sym(nwm1,s);
mW1_d = poly2sym(dwm1,s);
W1_ms = mW1_n/mW1_d;

% Calculation of W2(-s)
for k = 1:length(nw2)
    nwm2(k) = nw2(k)*(-1)^(length(nw2)-k);
end
for k = 1:length(dw2)
    dwm2(k) = dw2(k)*(-1)^(length(dw2)-k);
end
mW2_n = poly2sym(nwm2,s);
mW2_d = poly2sym(dwm2,s);
W2_ms = mW2_n/mW2_d;

den_E_gamma = W1_d*mW1_d*gamma_opt^2;
num_E_gamma = W1_n*mW1_n-den_E_gamma;
% E_gamma = minreal(sym2tf(num_E_gamma/den_E_gamma),EPSSS);
E_gamma = tf2sym(sym2tf(num_E_gamma/den_E_gamma));

a1=polyadd(conv(nw1,star(nw1)),-gamma_opt^2*conv(dw1,star(dw1)));
a1 = a1(find(a1,1):end);
b1=gamma_opt^2*conv(dw1,star(dw1)); 
b1 = b1(find(b1,1):end);

a2=polyadd(conv(nw2,star(nw2)),-gamma_opt^2*conv(dw2,star(dw2)));
a2 = a2(find(a2,1):end);
b2=gamma_opt^2*conv(dw2,star(dw2));
b2 = b2(find(b2,1):end);

A1=conv(a1,a2);
B1=conv(b1,b2);
A=polyadd(B1,-A1);
B=B1;
[dG,nG]=spec(A,B);
eta=roots(dw1);
n1=length(eta);
nF=conv(nG,poly(-eta));
dF=conv(dG,poly(eta));
hat_G_gamma = -poly2sym(nG,s)/poly2sym(dG,s)/gamma_opt;%%%!!!!!!
hat_G_gamma = minreal(sym2tf(hat_G_gamma),EPSSS);
hat_G_gamma = tf2sym(hat_G_gamma);
G_gamma = hat_G_gamma*W1_s;
G_gamma = minreal(sym2tf(G_gamma),EPSSS);
G_gamma = tf2sym(G_gamma);
num_F_gamma = poly2sym(nF,s);
den_F_gamma = poly2sym(dF,s);
% F_gamma = minreal(sym2tf(num_F_gamma/den_F_gamma),EPSSS);
F_gamma = tf2sym(sym2tf(num_F_gamma/den_F_gamma));
beta=roots(sym2poly(num_E_gamma));
beta_checker = beta(abs(real(beta))<=1e-5);
indices = [];
for k=1:length(beta_checker)
    if imag(beta_checker(k))==0
        indices = [indices k];
    end
end
beta_checker = beta_checker(indices(1:length(indices)/2));
beta=[beta((real(beta)>1e-5) | ((abs(real(beta))<=1e-5) & (imag(beta)>0))) beta_checker];
len_alpha = length(alpha);
if length(beta)==length(roots(dw1))
    M=[];
    for k=1:len_alpha,
       s = alpha(k);
       M=[M; alpha(k).^(0:n1+len_alpha-1) eval(Mn)*eval(F_gamma)*alpha(k).^[0:n1+len_alpha-1]];
    end;
    for k=1:n1,
       s = beta(k);
       M=[M; beta(k).^[0:n1+len_alpha-1] eval(Mn)*eval(F_gamma)*beta(k).^(0:n1+len_alpha-1)];
    end;
    for k=1:len_alpha,
       s = alpha(k);
       M=[M;(-alpha(k)).^[0:n1+len_alpha-1]*eval(Mn)*eval(F_gamma) (-alpha(k)).^[0:n1+len_alpha-1]];
    end;
    for k=1:n1,
       s = beta(k);
       M=[M; (-beta(k)).^[0:n1+len_alpha-1]*eval(Mn)*eval(F_gamma) (-beta(k)).^[0:n1+len_alpha-1]];
    end;
    
    [~,Suu,Vuu]=svd(M);
    LL=Vuu(:,length(M));
    Suu(length(M),length(M))
    if Suu(length(M),length(M)) < eps
        syms s;
        num_L_gamma=poly2sym(real(LL(2*(n1+len_alpha):-1:n1+len_alpha+1))+imag(LL(2*(n1+len_alpha):-1:n1+len_alpha+1)),s);
        den_L_gamma=poly2sym(real(LL(n1+len_alpha:-1:1))+imag(LL(n1+len_alpha:-1:1)),s);
        L_gamma = num_L_gamma/den_L_gamma;
        L_gamma = minreal(sym2tf(L_gamma),EPSSS);
        L_gamma = tf2sym(L_gamma);
        K_opt = 1/L_gamma;
        K_opt = minreal(sym2tf(K_opt),EPSSS);
        K_opt = tf2sym(K_opt);
        sym2zpk(K_opt)
        Ro = W1_d*W1_d/(num_E_gamma*Md);
        Ro = minreal(sym2tf(Ro),EPSSS);
        Ro = tf2sym(Ro);
        s=1e+45;
        dInf=real(eval(Ro)*eval(K_opt)/gamma_opt);
        syms s;
        M1 = mW1_d/W1_d;
        [Hn,Hd]=findHnHd(gamma_opt,K_opt,Ro,Mn,M1,hat_G_gamma,dInf);
        Co = 1/(1+Hn);
        lw=lw1:(lw2-lw1)/Lpoints:lw2;
        wval=10.^lw;
        
        try sym2tf(Hd)
            if (order_inputs.DO_APP_EVEN_FINITE==1)
                app_Hd = doApprox(Hd,wval,order_inputs.Hd_order,order_inputs.Hd_rel_order,Co/(1+Co*Hd));
                Hd_bode_title = ['Bode Plot of H_d(s) and Its ',num2str(order_inputs.Hd_order),' Order Approximation'];
            else
                app_Hd = Hd;
                Hd_bode_title = 'Bode Plot of H_d(s) Without Approximation (Finite Dimensional)';
            end
        catch
            app_Hd = doApprox(Hd,wval,order_inputs.Hd_order,order_inputs.Hd_rel_order,Co/(1+Co*Hd));
            Hd_bode_title = ['Bode Plot of H_d(s) and Its ',num2str(order_inputs.Hd_order),' Order Approximation'];
        end
        
        try sym2tf(No)
            if (order_inputs.DO_APP_EVEN_FINITE==1)
                app_No = doApprox(No,wval,order_inputs.No_order,order_inputs.No_rel_order,No);
                No_bode_title = ['Bode Plot of N_o(s) and Its ',num2str(order_inputs.No_order),' Order Approximation'];
            else
                app_No = No;
                No_bode_title = 'Bode Plot of N_o(s) Without Approximation (Finite Dimensional)';
            end
        catch
            app_No = doApprox(No,wval,order_inputs.No_order,order_inputs.No_rel_order,No);
            No_bode_title = ['Bode Plot of N_o(s) and Its ',num2str(order_inputs.No_order),' Order Approximation'];
        end
        
        term = (Co/(1+Co*app_Hd));
        Ca = 1/(dInf*gamma_opt^2)*term*hat_G_gamma/app_No;
        Copt = 1/(dInf*gamma_opt^2)*(Co/(1+Co*Hd))*hat_G_gamma/No;
%         Ca = 4.04*(1+2*0.37*s/0.8+s^2/0.8^2)/(s*(1+s/30)); % Article - K for h=0
        
        s = 1j*wval;
        ca = eval(Ca);
%         c=(1/(dInf*gamma_opt^2)).*(eval(Co)./(1+eval(Co).*eval(Hd))).*eval(hat_G_gamma)./eval(No);
        c = eval(Copt);
         
        psi_inf=sqrt(abs(eval(W1_s)).^2+abs(eval(W2_s).*eval(Plant).*c).^2)./abs(1+eval(Plant).*c);
        psi_inf2=sqrt(abs(eval(W1_s)).^2+abs(eval(W2_s).*eval(Plant).*ca).^2)./abs(1+eval(Plant).*ca);
        p=eval(Plant);
        S=(1+p.*c).^(-1);
        T=1-S;
        perMatrix=[eval(W1_s).*S; eval(W2_s).*T];
        per = sqrt(abs(perMatrix(1,:)).^2+abs(perMatrix(2,:)).^2);
        cmag=abs(c);
        cang = angle(c);
        hd_app_mag=abs(eval(app_Hd));
        hd_app_ang=angle(eval(app_Hd));
        hdmag=abs(eval(Hd));
        hdang=angle(eval(Hd));
        camag=abs(ca);
        caang = angle(ca);
        no_app_mag=abs(eval(app_No));
        no_app_ang=angle(eval(app_No));
        nomag=abs(eval(No));
        noang=angle(eval(No));
        pcr=real(p.*c);
        pci=imag(p.*c);
        pcr2=real(p.*ca);
        pci2=imag(p.*ca);
        figure(6);
        clf;
        subplot(2,1,1);
        semilogx(wval,20*log10(cmag));
        hold on;
        semilogx(wval,20*log10(camag),'r');
        grid on;
        ylabel('Magnitude (dB)');
        legend('C_o_p_t','C_a')
        title('Bode Plots of Controllers');
        set(gca,'FontSize',16);

        subplot(2,1,2);
        semilogx(wval,unwrap(cang)*180/pi);
        hold on;
        semilogx(wval,unwrap(caang)*180/pi,'r');
        grid on;
        xlabel('Frequency (rad/s)');
        ylabel('Phase (deg)');
        set(gca,'FontSize',16);
        set(gcf,'Name','Bode Plots of Controllers');
        
        
        figure(7);
        clf;
        subplot(1,1,1);
        plot(pcr,pci,pcr,-pci,'-.');
        hold on;
        arrow(pcr(1),pci(1),pcr(10),pci(10));
        plot(-1,0,'kx','LineWidth',4)
        grid on;
        title('Nyquist plot with C_{opt}(j\omega)');
        set(gca,'FontSize',16);
        set(gcf,'Name','Nyquist plot with Copt(jw)');
        
        %%%%%%
        figure(8);
        clf;
        subplot(1,1,1);
        plot(pcr2,pci2,pcr2,-pci2,'-.');
        hold on;
        arrow(pcr2(1),pci2(1),pcr2(10),pci2(10));
        plot(-1,0,'kx','LineWidth',4)
        grid on;
        title('Nyquist plot with C_a(j\omega)');
        set(gca,'FontSize',16);
        set(gcf,'Name','Nyquist plot with Ca(jw)');
        
        
        figure(9)
        clf;
        subplot(2,1,1);
        semilogx(wval,20*log10(hdmag));
        grid on; hold on;
        semilogx(wval,20*log10(hd_app_mag));
        ylabel('Magnitude (dB)');
        title(Hd_bode_title);
        legend('H_d','H_d_a');
        set(gca,'FontSize',16);
        
        subplot(2,1,2);
        semilogx(wval,unwrap(hdang)*180/pi);
        grid on;
        hold on;
        semilogx(wval,unwrap(hd_app_ang)*180/pi);
        xlabel('\omega (rad/s)');
        ylabel('Phase (deg)');
        set(gca,'FontSize',16);
        set(gcf,'Name','Bode Plot of Hd(s)');
        
        
        figure(10)
        clf;
        subplot(2,1,1);
        semilogx(wval,20*log10(nomag));
        grid on;
        hold on;
        semilogx(wval,20*log10(no_app_mag));
        ylabel('Magnitude (dB)');
        title(No_bode_title);
        legend('N_o','N_{oa}');
        set(gca,'FontSize',16);
        
        subplot(2,1,2);
        semilogx(wval,unwrap(noang)*180/pi);
        grid on;
        hold on;
        semilogx(wval,unwrap(no_app_ang)*180/pi);
        xlabel('\omega (rad/s)');
        ylabel('Phase (deg)');
        set(gca,'FontSize',16);
        set(gcf,'Name','Bode Plot of No(s)');
        
        figure(11);
        clf;
        semilogx(wval,(psi_inf),'b');
        hold on
        semilogx(wval,(psi_inf2),'r');
        gamma_app = max(psi_inf2)
        grid;
        xlabel('Frequency (rad/s)');
        ylabel('Performance level');
        set(gca,'fontsize', 22)
        title('Performance plot');
        legend('C_o_p_t','C_a');
        set(gca,'FontSize',16);
        set(gcf,'Name','Performance plot');
        
        warning = 'Correct GAMMA_OPT value';
    else
        warning = 'Incorrect GAMMA_OPT value'
    end
else
	warning = 'Non-generic case and/or incorrect GAMMA_OPT value'
end