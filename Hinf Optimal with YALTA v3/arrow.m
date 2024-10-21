function arrow(x1,y1,x2,y2)
z0=x2+i*y2;
z=(x2-x1)+i*(y2-y1);
ax=axis;
z=z/abs(z)*max(ax(2)-ax(1),ax(4)-ax(3))/40;
z1=z0-z*exp(i*pi/6);
z2=z0-z*exp(-i*pi/6);
plot([real(z1),x2],[imag(z1),y2]);
plot([real(z2),x2],[imag(z2),y2]);

