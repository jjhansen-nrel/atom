function f=exp23B_12(x1,y1,x2,y2,l);
y=y1-y2;
x=x1-x2;
r=sqrt(x.^2+y.^2);
ind0=(r==0);
ind=(~ind0);
f(ind0)=0;
r1=r(ind);
f(ind)=exp(-(r1/l).^(2/3)).*(x(ind).*y(ind)./(3*r1.^(4/3)*l^(2/3)));