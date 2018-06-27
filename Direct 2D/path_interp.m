function s_out = path_interp(s_in, xbeg, ybeg, azi, sgrid, xaxis, yaxis)
%PATH_INT  Interpolates a quantity along a path through a 2D domain.

x = xbeg + s_in*cos(azi);
y = ybeg + s_in*sin(azi);
s_out = interp2(xaxis, yaxis, sgrid, x, y);
