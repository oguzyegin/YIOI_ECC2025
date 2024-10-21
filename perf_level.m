function [psi,g] = perf_level(W1,W2,Pdf,Ms,C,w)
W1 = tf2sym(tf(W1));
W2 = tf2sym(tf(W2));
Pdf = tf2sym(tf(Pdf));
C = tf2sym(tf(C));

s = 1j*w;

S = 1./(1+eval(Pdf).*eval(Ms).*eval(C));
T = eval(Pdf).*eval(Ms).*eval(C)./(1+eval(Pdf).*eval(Ms).*eval(C));

psi = sqrt(abs(eval(W1).*S).^2+abs(eval(W2).*T).^2);
g = max(psi);