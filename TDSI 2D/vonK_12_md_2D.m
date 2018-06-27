function f=vonK_12_md_2D(x_vec,path,Kv,xi,yj);
x=path.x0+x_vec*path.cos;
y=path.y0+x_vec*path.sin;
f=vonKB_12(x,y,xi,yj,Kv);