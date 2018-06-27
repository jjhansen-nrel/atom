function f=expB_12(x1,y1,x2,y2,l);
y=y1-y2;
x=x1-x2;
r=sqrt(x.^2+y.^2);
ind0=(r==0);
ind=(~ind0);
f(ind0)=0;
f(ind)=exp(-r(ind)/l).*(x(ind).*y(ind)./(2*r(ind)*l));