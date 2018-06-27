function f=exp23B_22(x1,y1,x2,y2,l);
x=(x1-x2).^2;
r=sqrt((y1-y2).^2+x);
ind0=(r==0);
i=(~ind0);
f(ind0)=1;
r1=r(i);
f(i)=exp(-(r1/l).^(2/3)).*(1-x(i)./(3*r1.^(4/3)*l^(2/3)));