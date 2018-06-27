function f=expB_22(x1,y1,x2,y2,l);
x=(x1-x2).^2;
r=sqrt((y1-y2).^2+x);
ind0=(r==0);
i=(~ind0);
f(ind0)=1;
f(i)=exp(-r(i)/l).*(1-x(i)./(2*r(i)*l));