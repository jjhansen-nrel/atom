function f=expB_11(x1,y1,x2,y2,l);
y=(y1-y2).^2;
r=sqrt((x1-x2).^2+y);
ind=(r==0);
f(ind)=1;
r1=r(~ind);
f(~ind)=exp(-r1/l).*(1-y(~ind)./(2*r1*l));