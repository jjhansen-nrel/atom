function f=exp23_12_md_2D(x_vec,path,lv,xi,yj);
x=path.x0+x_vec*path.cos;
y=path.y0+x_vec*path.sin;
f=exp23B_12(x,y,xi,yj,lv);