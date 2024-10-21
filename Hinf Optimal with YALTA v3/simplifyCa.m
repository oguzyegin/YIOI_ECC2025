function Ca = simplifyCa(Ca)
syms s
newNum = s/s;
newDen = s/s;
[num,den]=numden(Ca);
num = sym2poly(num);
a0=num(1);
den = sym2poly(den);
b0=den(1);
rootOfNum = roots(num);
rootOfDen = roots(den);
sameRoots=[]; indices=[];
count=0;
for k=1:length(rootOfNum)
    currentRoot=rootOfNum(k);
    sameRoot=rootOfDen(abs(rootOfDen-currentRoot)<1e-8);
    if isempty(sameRoot)
        count=count+1;
        rootsOfNum(count)=currentRoot;
    end
end
for k=1:count
    newNum = newNum*(s-rootsOfNum(k)); 
end
count=0;
for k=1:length(rootOfDen)
    currentRoot=rootOfDen(k);
    sameRoot=rootOfDen(abs(rootOfNum-currentRoot)<1e-8);
    if isempty(sameRoot)
        count=count+1;
        rootsOfDen(count)=currentRoot;
    end
end
for k=1:count
    newDen = newDen*(s-rootsOfDen(k)); 
end
newNum=sym2tf(newNum);
newDen=sym2tf(newDen);
Ca = (a0/b0)*newNum/newDen;