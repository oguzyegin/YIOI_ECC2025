function s=addpoly(p,q)
np=length(p);
nq=length(q);
n=max(np,nq);
if n==np
    s=[zeros(1,n-nq) q]+p;
else
    s=[zeros(1,n-np) p]+q;
end



