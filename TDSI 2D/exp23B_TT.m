function f=exp23B_TT(x1,y1,x2,y2,l);
r=sqrt((x1-x2).^2+(y1-y2).^2);
f=exp(-(r/l).^(2/3));