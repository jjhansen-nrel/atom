function y=I3(my,e1,e2,e3,sin_phi,cos_phi,y2,a,b,c,Li,lv,sigmav2);

% my=path1.y0-y2-a*path1.sin;
% e1=exp(-(Li+a)/lv.*(Li+a)/lv);
% e2=exp(-a.*a/lv/lv);
% e3=exp(-b.*b/(lv*lv));

y1=sin_phi*cos_phi*sigmav2*lv*e3.* ...
    (-(Li+a)/(2*lv).*e1+0.25*sqrt(pi)*erf((Li+a)/lv)+...
    a/(2*lv).*e2-0.25*sqrt(pi)*erf(a/lv));

y2=0.5*sigmav2*(cos_phi*cos_phi-sin_phi*sin_phi)*e3.*b.*(e1-e2);

y3=c/(lv*lv).*I1(a,b,Li,lv,sigmav2);

y=y1+y2+y3;