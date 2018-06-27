function f=exp_T_md_2D(x_vec,path,lT,xi,yj);
x=path.x0+x_vec*path.cos;
y=path.y0+x_vec*path.sin;
f=expB_TT(x,y,xi,yj,lT);