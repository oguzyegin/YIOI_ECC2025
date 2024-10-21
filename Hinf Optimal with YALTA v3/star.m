function q=star(p)

for k=1:length(p),
    q(k)=p(k)*(-1)^(length(p)+k);
end;
