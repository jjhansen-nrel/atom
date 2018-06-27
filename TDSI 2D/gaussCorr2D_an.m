function f=gaussCorr2D_an(x_vec,path1,path2,c0,T0,lT,lv,sigmaT2,sigmavx2,sigmavy2)
% Elements of R_dd matrix for Gaussian covariance functions. This function
% is used in the numerical integration in the calculation of R_dd elements;

x2=path2.x0+x_vec*path2.cos;
y2=path2.y0+x_vec*path2.sin;

Li=path1.length;
a=path1.cos*(path1.x0-x2)+path1.sin*(path1.y0-y2);
b=path1.sin*(path1.x0-x2)-path1.cos*(path1.y0-y2);
my=path1.y0-y2-a*path1.sin;
mx=path1.x0-x2-a*path1.cos;
e1=exp(-(Li+a)/lv.*(Li+a)/lv);
e2=exp(-a.*a/lv/lv);
e3=exp(-b.*b/(lv*lv));
dee3=e3.*(e1-e2);
bracee3=e3.*(-(Li+a)/(2*lv).*e1+0.25*sqrt(pi)*erf((Li+a)/lv)+...
    a/(2*lv).*e2-0.25*sqrt(pi)*erf(a/lv));
sin_phi=path1.sin;
cos_phi=path1.cos;
c=a.*b*(cos_phi*cos_phi-sin_phi*sin_phi)-a.*a*sin_phi*cos_phi+(path1.x0-x2).*(path1.y0-y2);

sigmav2=sqrt(sigmavx2*sigmavy2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% I2
J1=I1(a,b,Li,lv,sigmavx2);
J22=-my.*my/(lv*lv).*J1;
J21=-sin_phi*sin_phi*sigmavx2*lv*bracee3;
J23=sin_phi*sigmavx2*dee3.*my;
J2=J21+J22+J23;
I2=J1+J2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I3
y1=sin_phi*cos_phi*sigmav2*lv*bracee3;
y2=0.5*sigmav2*(cos_phi*cos_phi-sin_phi*sin_phi)*dee3.*b;
y3=c/(lv*lv).*I1(a,b,Li,lv,sigmav2);
I3=y1+y2+y3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I4
z1=I1(a,b,Li,lv,sigmavy2);
z22=-mx.*mx/(lv*lv).*z1;
z21=-cos_phi*cos_phi*sigmavy2*lv*bracee3;
z23=cos_phi*sigmavy2*dee3.*mx;
z2=z21+z22+z23;
I4=z1+z2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(c0)==1
    c0(1:2)=c0;
end
if length(T0)==1
    T0(1:2)=T0;
end
f=c0(1)*c0(2)/(4*T0(1)*T0(2))*I1(a,b,Li,lT,sigmaT2)+path1.cos*path2.cos*I2+...
    (path1.cos*path2.sin+path1.sin*path2.cos)*I3+...
    path1.sin*path2.sin*I4;
