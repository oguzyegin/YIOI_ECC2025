function q=pspec(p)

EPS=1e-12;

z=roots(p);
zn=z(find(real(z)<-EPS));
zi=z(find((abs(real(z))) <= EPS));

zx=[];
while (length(zi) > 0),
      zx=[zx; zi(1)];
      lz=length(zi);
      [x,n]=min(abs(zi(2:lz)-zi(1)*ones(lz-1,1)));
      n=n+1;
      zi=zi([2:n-1,n+1:lz]);
end;

q=poly([zn; zx]);
