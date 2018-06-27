function y=I2(my,e1,e2,e3,sin_phi,y2,a,b,Li,lv,sigmavx2);

% my=path1.y0-y2-a*path1.sin;
% e1=exp(-(Li+a)/lv.*(Li+a)/lv);
% e2=exp(-a.*a/lv/lv);
% e3=exp(-b.*b/(lv*lv));

J1=I1(a,b,Li,lv,sigmavx2);

J22=-my.*my/(lv*lv).*J1;

J21=-sin_phi*sin_phi*sigmavx2*lv*e3.* ...
    (-(Li+a)/(2*lv).*e1+0.25*sqrt(pi)*erf((Li+a)/lv)+...
    a/(2*lv).*e2-0.25*sqrt(pi)*erf(a/lv));

J23=sin+phi*sigmavx2*e3.*my.*(e1-e2);

J2=J21+J22+J23;

y=J1+J2;