function [newFunc,foundFunc]=pspec_sFunction(sFunc)   
syms s;
polyFunc = sym2poly(sFunc);
a0 = polyFunc(1);
polyFunc = polyFunc/a0;
if length(polyFunc)==1
	newFunc=(polyFunc);
else
	allRoots = roots(polyFunc);
    LHP = allRoots(real(allRoots)<=0);
    for k=1:length(LHP);
        if k==1
            newFunc = s-LHP(k);
            otherFunc = s+LHP(k);
        else
            newFunc = newFunc*(s-LHP(k));
            otherFunc = otherFunc*(s+LHP(k));
        end
    end
    newFunc = sqrt(a0)*newFunc;
    foundFunc = sqrt(a0)*newFunc*otherFunc;
end