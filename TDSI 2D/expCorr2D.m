function f=expCorr2D(x_vec,y_scal,path1,path2,c0,T0,lT,lv,sigmaT2,sigmavx2,sigmavy2)

x1=path1.x0+x_vec*path1.cos;
y1=path1.y0+x_vec*path1.sin;
x2=path2.x0+y_scal*path2.cos;
y2=path2.y0+y_scal*path2.sin;
sigmav2=sqrt(sigmavx2*sigmavy2);
f=c0(1)*c0(2)/(4*T0(1)*T0(2))*sigmaT2*expB_TT(x1,y1,x2,y2,lT)+path1.cos*path2.cos*sigmavx2*expB_11(x1,y1,x2,y2,lv)+...
    (path1.cos*path2.sin+path1.sin*path2.cos)*sigmav2*expB_12(x1,y1,x2,y2,lv)+...
    path1.sin*path2.sin*sigmavy2*expB_22(x1,y1,x2,y2,lv);
