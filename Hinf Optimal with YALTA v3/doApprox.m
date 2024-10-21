% % Returns approximation of TF by using given order, relative order and weight
function result = doApprox(mySym,wval,approx_order,rel_order,weight_sym)
s = 1j*wval;
if nargin<5
    WFd = 1;
else
    WFd = frd(eval(weight_sym),wval);
end
data = eval(mySym);
dc_gain = abs(data(1));
data = data/dc_gain;
FFd=frd(data,wval);
if isempty(rel_order)
    Ff=fitfrd(FFd,approx_order,[],WFd);
%     Ff=tfest(FFd,approx_order,approx_order);
else
    Ff=fitfrd(FFd,approx_order,rel_order,WFd);
%     Ff=tfest(FFd,approx_order,approx_order-rel_order);
end
% Ff = tf(Ff);
% result_n = Ff.num{1}; result_d = Ff.den{1};
[result_n,result_d] = ss2tf(Ff.a,Ff.b,Ff.c,Ff.d);
syms s;
result = dc_gain*poly2sym(result_n,s)/poly2sym(result_d,s);