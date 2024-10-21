% Used https://nl.mathworks.com/help/symbolic/inv.html
clc
clear all

n=6;

% BEST ORD 6 h=.75 k=2 Gamma_opt=10.518
% % moments
% sigma1=2.375;sigma2=0.685;
% f1 = 7; f2 =12.439;
% s1=0; s2=4.3; 
% s3=sigma1-f1*1i; s4=sigma1+f1*1i; 
% s5=sigma2-f2*1i; s6=sigma2+f2*1i;
%Gamma_opt=10.518 
%s1=0;s2=4.3;s34=2.375+/-7j;s56=4.3+/-12.439j
% Kr6_tf =
% 
%      0.5732 s^5 + 137.8 s^4 - 384.1 s^3 + 5.263e04 s^2 - 1.039e05 s + 3.044e06
%   --------------------------------------------------------------------------------
%   s^6 + 33.7 s^5 + 750.2 s^4 + 1.042e04 s^3 + 9.543e04 s^2 + 5.248e05 s + 1.335e06
%poles  -3.8381 +     14.276i
      % -3.8381 -     14.276i
      % -7.0695 +     2.7203i
      % -7.0695 -     2.7203i
      % -5.9436 +     8.4326i
      % -5.9436 -     8.4326i

% moments
sigma1=2.375;sigma2=0.685;
f1 = 7; f2 =12.439;
s1=0; s2=4.3; 
s3=sigma1-f1*1i; s4=sigma1+f1*1i; 
s5=sigma2-f2*1i; s6=sigma2+f2*1i;

% moments
%f1 = 10; f2 =18;
%s1=0; s2=12; s3=f1*1i; s4=-f1*1i; s5=f2*1i; s6=-f2*1i;

% Define complex variable s
s = sym('s'); % Declare s as a symbolic variable
sx = sym('sx'); % Declare sx as a symbolic variable


% Impose that the G variables are real
%assume([g11, g12, g13, g14, g15, g16], 'real');

syms g11 g12 g13 g14 g15 g16 real

%% Controller and  Hopt transfer functions
%sx=tf('s');
% this is for a=1
aa=1;
hhat=0.75; khat=2;
phopt=-0.081637; gopt=8.96475; 
kghat=khat*sqrt((2/khat)-gopt^(-2));

kghat=khat*sqrt((2/khat)-gopt^(-2));

% Define CpiTF, HoptTF, CoptTF and GoptTF
CpiTF=gopt*(s-phopt)/(khat*s);

HoptTF(s)= (s*(phopt*s-gopt^(-2))+(s^2+gopt^(-2))+...
      (1/khat)*( (s+phopt)*(kghat*s+1)+gopt*s*(phopt-s)*exp(-hhat*s)))/...
           ( (s^2+gopt^(-2))*(s-1) );
CoptTF=CpiTF/(1+HoptTF);
GoptTF=CoptTF*exp(-hhat*s)/(s-1);
Sopt=1/(1+GoptTF);
Topt=1-Sopt;
ClTF=sqrt(abs(Sopt/s)^2 + abs(khat*s*Topt)^2);% optimal is flat
        % the value is constant at gopt=4.944 as expected...
        %W1SMM(k)=abs(SMM(k)/ss);
        %W2TMM(k)=abs(khat*s*TMM(k));

% Define K(s) which is HoptTF
% Compute the derivative of K with respect to s
K(s) =  (s*(phopt*s-gopt^(-2))+(s^2+gopt^(-2))+...
      (1/khat)*( (s+phopt)*(kghat*s+1)+gopt*s*(phopt-s)*exp(-hhat*s)))/...
           ( (s^2+gopt^(-2))*(s-1) );
dK(s) = diff(K, s);

%%% Calculate moments of order zero and one
% Substitute moments into K(s) and the derivative dK(s)
% note that one needs to use double command
%%% Calculate moments of order zero and one
% Convert the symbolic substitutions and ensure they are double precision
eta = double(subs(K(s), s, [s1, s2, s3, s4, s5, s6]));
gamma = double(subs(dK(s), s, [s1, s2, s3, s4, s5, s6]));

disp('We are matching 6 derivatives.');
disp('The interpolation points (col 1) the moments of ord 0 and 1 (col 2 and 3) are:');
disp([s1 eta(1) gamma(1); s2 eta(2) gamma(2); s3 eta(3) gamma(3); s4 eta(4) gamma(4); s5 eta(5) gamma(5); s6 eta(6) gamma(6)]);


% Define the matrix S and L
S = [s1, 0, 0, 0, 0, 0; 
     0, s2, 0, 0, 0, 0; 
     0, 0, sigma1, f1, 0, 0;
     0, 0, -f1, sigma1, 0, 0;
     0, 0, 0, 0, sigma2, f2;
     0, 0, 0, 0, -f2, sigma2]
L = [1, 1, 0, sqrt(2), 0, sqrt(2)]

%Define G and F, where F = S - G*L
G = [g11; g12; g13; g14; g15; g16];
F = S - G * L

% Define the matrix H as a 1x6 row vector using the real and imaginary parts
%H = [eta0, eta1, -imag(eta3), real(eta3), -imag(eta5), real(eta5)]
H = [eta(1), eta(2), imag(eta(3))*sqrt(2), real(eta(3)*sqrt(2)), imag(eta(5))*sqrt(2), real(eta(5))*sqrt(2)];
Hg=double(H); % Let Hg be of type double
disp('H and Hg are the same as it should be');
disp(H);
disp(Hg);

% Compute Krg as (H * inv(s - F) * G)
inv_sF = inv(s * eye(6) - F);
Krg = H * inv_sF * G;

% Extract the (1,1) element of the result
% Display the result
Krg_11 = Krg(1, 1);
disp('The (H * inv(s - F) * G)(1,1) is:');
disp(simplify(Krg_11));
dKrg_11(s) = diff(Krg_11, s);


% Compute derivatives and evaluate at specific points
% Display the results
m = subs(dKrg_11(s), s, [s1, s2, s3, s4, s5, s6]);
disp('The values of m are:');
disp(m);


% Define the equations based on m0, m1, m2, m3, m4, and m5
% eq1 = double(m0) - gamma0; % Equation for m0
eqs = m - gamma;

% Create a function for fsolve
% Convert the equations into a function handle for fsolve using matlabFunction
fsolve_func = matlabFunction(eqs, 'Vars', {[g11, g12, g13, g14, g15, g16]});

% Initial guess for the variables [g11, g12, g13, g14, g15, g16]
initial_guess = [0.5 ,0.5 ,0.5 , 0.5 ,0.5 ,0.5]; % Change this if you have a better guess

% Set options with tolerances for controlling the precision
options = optimoptions('fsolve', 'TolFun', 1e-12, 'TolX', 1e-12, 'Display', 'iter');

% Use fsolve to solve the equations
[gg, fval, exitflag] = fsolve(fsolve_func, initial_guess, options);

% Display the results
Gg=double(real(gg'));
disp('The solution for g11, g12, g13, g14, g15, g16 is:');
disp(Gg);

% Find Fg
F_subst = double(S - Gg.* L) % Substitute the values in Gg into F
Fg_imaginar=imag(F_subst) % verify Fg is zero 
Fg=real(F_subst) % to aavoit unnecessary errors due to zero imag part
disp('The eigenvalues of Fg are:');
eig(Fg)

%% Create tranffer function for Kr6
sys_ss = ss(Fg, Gg, H, 0);% Pas 1. Create the state-space model
Kr6_tf = tf(sys_ss);% Pas 2. Convert the state-space model to a transfer function
disp('The transfer function is:');% Pas 3. Display the transfer function
Kr6_tf

% Extract numerator and denominator
[num, den] = tfdata(Kr6_tf, 'v'); % 'v' option gives the result as vectors

% % Display numerator and denominator
% disp('Numerator of Kr6_tf:');
% disp(num);
% disp('Denominator of Kr6_tf:');
% disp(den);

% Construct the transfer function in symbolic form
Kr6_symbolic = poly2sym(num, s) / poly2sym(den, s);

% Display the symbolic transfer function
disp('The symbolic transfer function is:');
disp(Kr6_symbolic);


%%%%% De calculat aproximantul Kr6 pt evaluari in momente s1 etc
%% Compute Krg as (H * inv(s - F) * G)
Kr6_in_s(s) = simplify(H * (inv(s * eye(6) - Fg)) * Gg)


%De pus restul check the moments
HoptTF_at_s1 = double(subs(HoptTF(s), s, s1));
Kr6_in_s_at_s1 = subs(Kr6_in_s(s), s, s1);
HoptTF_at_s2 = double(subs(HoptTF(s), s, s2));
Kr6_in_s_at_s2 = subs(Kr6_in_s(s), s, s2);
HoptTF_at_s3 = double(subs(HoptTF(s), s, s3));
Kr6_in_s_at_s3 = double(subs(Kr6_in_s(s), s, s3));
HoptTF_at_s4 = double(subs(HoptTF(s), s, s4));
Kr6_in_s_at_s4 = double(subs(Kr6_in_s(s), s, s4));
HoptTF_at_s5 = double(subs(HoptTF(s), s, s5));
Kr6_in_s_at_s5 = double(subs(Kr6_in_s(s), s, s5));
HoptTF_at_s6 = double(subs(HoptTF(s), s, s6));
Kr6_in_s_at_s6 = double(subs(Kr6_in_s(s), s, s6));


% Display the result
disp('The value of HoptTF(s1)-Kr6_symbolic_at_s1 is:');
disp(HoptTF_at_s1-Kr6_in_s_at_s1);
disp(HoptTF_at_s2-Kr6_in_s_at_s2);
disp(HoptTF_at_s3-Kr6_in_s_at_s3);
disp(HoptTF_at_s4-Kr6_in_s_at_s4);
disp(HoptTF_at_s5-Kr6_in_s_at_s5);
disp(HoptTF_at_s6-Kr6_in_s_at_s6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Bode plots HKr6 and Hopt
% Frequency range
om = logspace(-1, 3, 2000);
s_freq = 1i * om; % Define s = jω for the frequency response

% Convert K(s) to a transfer function for evaluation Hopt
K_tf = matlabFunction(K); % Convert symbolic K(s) to a function handle

% Evaluate K(s) over the frequency range      
K_response = K_tf(s_freq);       % Evaluate K(s) at the frequency points
mag_K = abs(K_response);    % Magnitude of K(s)
phase_K = angle(K_response) * (180/pi); % Phase of K(s) in degrees
Hopt = K_response;


% Get the complex frequency response of Kr6_tf at the defined frequencies
[magKr6_tf, phaseKr6_tf, w] = bode(Kr6_tf, om);

% Bode function returns a 3D array; we need to squeeze the magnitude array
magKr6_tf = squeeze(magKr6_tf); % Extract the magnitude
phaseKr6_tf = squeeze(phaseKr6_tf); % Extract the phase


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Evaluate CpiTF, CMM, and ClTFMM for the obtained approximant
Kr6_response = matlabFunction(Kr6_symbolic);
HKrMM = Kr6_response(s_freq);% Assign HKrMM as Kr6_response

CpiTF_response = matlabFunction(CpiTF); % Convert symbolic CpiTF(s) to a function handle
Cpi=CpiTF_response(s_freq); % Evaluate CpiTF(s) at the frequency points using the function handle
% Define CMM as Cpi / (1 + HKrMM)
CMM = Cpi ./ (1 + HKrMM); % Element-wise division for each frequency point
GMM = CMM .* exp(-hhat * s_freq) ./ (s_freq - 1);
SMM = 1./(1+GMM);
TMM = 1-SMM;
ClTFMM=sqrt(abs(SMM ./s_freq).^2 + abs(khat * s_freq.*TMM).^2);

%Evaluate ClTF over the frequency range
CoptTF_resp=Cpi ./(1 + Hopt);
GoptTF_resp=CoptTF_resp .* exp(-hhat*s_freq) ./ (s_freq-1);
Sopt_resp=1 ./ (1+GoptTF_resp);
Topt_resp=1-Sopt_resp;
ClTFs_resp=sqrt(abs(Sopt_resp ./ s_freq).^2 + abs(khat*s_freq .* Topt_resp).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%0.75 k=2 Ord=6
% Pade from report; 
% HPade6 = -2.5703*(s - 131.5)*(s^2 - 2.271*s + 71.82)*(s^2 - 4.61*s + 351.3)/...
% (s^2 + 22.66*s + 133.7)*(s^2 + 19.92*s + 148.3)*(s^2 + 13.42*s + 188.5);


% s1=?, s2=?, s34=+-?j, s56=+-?j 
%Systune% Ord 6 cu gain 
% HSysTune6 = (0.044446*(s^2-1.962*s + 71.45)*(s^2-1.163*s + 289.1)*(s^2 + 22.11*s + 1293))/...
%     ((s^2 + 4.943*s + 15.71)*(s^2 + 2.768*s + 107.5)*(s^2 + 1.289*s + 321.7));

%systune from report
HSysTune6 = 1.1429*(s + 34.63)*(s^2 - 2.542*s + 78.55)*(s^2 - 2.264*s + 286.7)/...
((s^2 + 4.858*s + 14.78)*(s^2 + 2.857*s + 98.74)*(s^2 + 1.117*s + 291.6));

HKrST_response = matlabFunction(HSysTune6);
   HKrST=HKrST_response(s_freq);

    CST = Cpi./(1+HKrST);
    GST = CST.*exp(-hhat * s_freq) ./ (s_freq-1);
    SST = 1./(1+GST);
    TST = 1-SST;
    ClTFST = sqrt(abs(SST./s_freq).^2+abs(khat*s_freq.*TST).^2);   

% IRKA from report; 
 HIRKA6 = 0.61126*(s + 166.9)*(s^2 - 2.268*s + 71.97)*(s^2 - 2.298*s + 282.5)/...
 ((s^2 + 9.055*s + 30.24)*(s^2 + 7.423*s + 109)*(s^2 + 5.272*s + 277.6));

 HIRKA6_response = matlabFunction(HIRKA6);
   HIRKA=HIRKA6_response(s_freq);

    CIRKA = Cpi./(1+HIRKA);
    GIRKA = CIRKA.*exp(-hhat * s_freq) ./ (s_freq-1);
    SIRKA = 1./(1+GIRKA);
    TIRKA = 1-SIRKA;
    ClTFIRKA = sqrt(abs(SIRKA./s_freq).^2+abs(khat*s_freq.*TIRKA).^2);

 
%%%%%%%%%%%%%%%%%%%%%%
%PLOTS to compare 3
%%%%%%%%%%%%%%%%%%%%%%
% Magnitude plot
% Create a figure with two subplots
figure(1);

% Plot the absolute value / the magnitude
subplot(2, 1, 1);
semilogx(w, mag_K,'b',w, magKr6_tf,'r',w, abs(HKrST),'k-.',w, abs(HIRKA),'g--');
xlabel('Frequency (rad/s)');
ylabel('|TF(j\omega)|');
title('Magnitude of MM6_{tf}, SysTune_{tf}, IRKA_{tf} and HKopt_{tf}');
grid on;

% Plot the phase
subplot(2, 1, 2);
semilogx(w,phase_K,'b',w, phaseKr6_tf, 'r',w, angle(HKrST),'k-.',w, angle(HIRKA),'g--');
xlabel('Frequency (rad/s)');
ylabel('Phase (degrees)');
grid on;


% %% Error plot: Hopt-HoptApprox  abs(HSysTune6-Hopt) and abs(HIRKA-Hopt)
figure(11)
%loglog(om,abs(HKrMM-Hopt),'r')
loglog(om,abs(HKrMM-Hopt),'r',om,abs(HKrST-Hopt),'k-.',om,abs(HIRKA-Hopt),'g--')
title({'h=.75, k=2, \nu=6: Error plot ', '|MM6-Hopt|-red'},' |SysTune6-Hopt|-black, |IRKA6-Hopt|-green')

% %% gamma plot ClTFMM    ClTFST and ClTFIRKA
 figure(12) % gamma
 loglog(om,ClTFs_resp,'b',om,ClTFMM,'r',om,ClTFST,'k-.',om,ClTFIRKA,'g--')
% %semilogx(om,ClTFMM,'r')
% semilogx(om,ClTF,'b',om,ClTFMM,'r',om,ClTFST,'k-.')
 title({'h=.75, k=2, \nu=6:  \gamma plot','Copt-blue (peak=8.96), MM6-red (peak=10.518)'}, 'SysTune66-black (peak=11.0852), IRKA6-black (peak=11.0166)')

% figure(13)
% loglog(om,abs(ClTFMM-ClTF),'r',om,abs(ClTFST-ClTF),'k-.')
% title('Error plot (h=.75 k=2): ClTFST-red ord 6')

%%%%%%%%%%%%%%%%%%%%%%
%Figures to comapre only Opt MM6 and SysTune6
%%%%%%%%%%%%%%%%%%%%%%
% % Magnitude plot
% % Create a figure with two subplots
% figure(1);
% 
% % Plot the absolute value / the magnitude
% subplot(2, 1, 1);
% semilogx(w, mag_K,'b',w, magKr6_tf,'r',w, abs(HKrST),'k-.');
% xlabel('Frequency (rad/s)');
% ylabel('|TF(j\omega)|');
% title('Magnitude of MM6_{tf} SysTune_{tf} and HKopt_{tf}');
% grid on;
% 
% % Plot the phase
% subplot(2, 1, 2);
% semilogx(w,phase_K,'b',w, phaseKr6_tf, 'r',w, angle(HKrST),'k-.');
% xlabel('Frequency (rad/s)');
% ylabel('Phase (degrees)');
% grid on;
% 
% 
% % %% Error plot: Hopt-HoptApprox and abs(HSysTune6-Hopt)
% figure(11)
% %loglog(om,abs(HKrMM-Hopt),'r')
% loglog(om,abs(HKrMM-Hopt),'r',om,abs(HKrST-Hopt),'k-.')
% title('Error plot (h=.75 k=2): HKrMM-red ord 6, ClTFST6-black')
% 
% % %% gamma plot ClTFMM  To Be Added ClTFST
%  figure(12) % gamma
%  loglog(om,ClTFs_resp,'b',om,ClTFMM,'r',om,ClTFST,'k-.')
% % %semilogx(om,ClTFMM,'r')
% % semilogx(om,ClTF,'b',om,ClTFMM,'r',om,ClTFST,'k-.')
%  title({'Gamma (h=.75 k=2):',' Copt-blue (peak=8.96)', 'ord 6 ClTFMM-red (peak=10.518)', 'ClTFST6-black (peak=10.852)'})
% 
% % figure(13)
% % loglog(om,abs(ClTFMM-ClTF),'r',om,abs(ClTFST-ClTF),'k-.')
% % title('Error plot (h=.75 k=2): ClTFST-red ord 6')