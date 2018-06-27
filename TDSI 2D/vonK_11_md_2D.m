function f=vonK_11_md_2D(x_vec,path,K0,xi,yj);
x=path.x0+x_vec*path.cos;
y=path.y0+x_vec*path.sin;
f=vonKB_11(x,y,xi,yj,K0);