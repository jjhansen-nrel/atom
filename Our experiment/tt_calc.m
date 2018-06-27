function [tt,dtt,Li]=tt_calc(TT,VX,VY,xy,xv,yv,sigman,S,R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC=343*sqrt(TT/293);
xsrc=xy(1:S,1);
ysrc=xy(1:S,2);
xrcv=xy(S+1:end,1);
yrcv=xy(S+1:end,2);
Ntime=size(TT,3);
c0=zeros(Ntime,1);
dtt=zeros(S*R,Ntime);
tt=zeros(S*R,Ntime);
Li=zeros(S*R,1);

  h_bar=waitbar(0,'Travel times calculation','name','0%');
for i=1:Ntime
  C=CC(:,:,i);
  c0(i)=mean(C(:));
  vx=VX(:,:,i);
  vy=VY(:,:,i);
  tt_pert = ForwardProb(xsrc, ysrc, xrcv, yrcv, xv, yv, C', vx', vy', c0(i), sigman);
  tt_pert1=tt_pert';
  dtt(:,i)=tt_pert1(:);
  waitbar(i/Ntime,h_bar);
  set(h_bar,'name',sprintf('%d%%',round(i/Ntime*100)));
end
close(h_bar);
for i=1:S % loop over all sources
    for j=1:R % loop over all receivers
        b=xy(S+j,1)-xy(i,1); % difference of x between the j-th receiver and i-th source 
        a=xy(S+j,2)-xy(i,2); % difference of y between the j-th receiver and i-th source
        Li(R*(i-1)+j,1)=sqrt(a*a+b*b);
    end
end
for t=1:Ntime
    tt(:,t)=Li/c0(t)+dtt(:,t);
end
