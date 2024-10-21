function [string]=nicePrint(Param)
syms s;
try
    Param = simplifySym(Param);
    [num,den] = numden(Param);
    num = sym2poly(num);
    den = sym2poly(den);
    [r,p,k]=residue(num,den);
    string='';
    if k~=0
        fin=length(r)+1;
    else
        fin=length(r);
    end
    lastEqn=0*s;
    for n=1:fin
        if n~=1
            if n>length(r)
                term=k*s/s;
            else
                term = r(n)/(s-p(n));
            end
        else
            term = r(n)/(s-p(n));
        end
        lastEqn = lastEqn+term;
    end
    str=evalc('pretty(vpa(lastEqn,4))');
    string=[string,str];
catch
    Param = simplify(Param);
    string=evalc('pretty(vpa(Param,2))');
end