function [meanField_T,meanField_V]=spatialMean_r(S,R,xy,tt,Lo,Hi,showFig)
% spatialMean_r does the same job for reciprocal transmission arrays as
% spatialMean does for non-reciprocal ones (see description of
% spatialMean). The difference is that spatialMean_r separates the travel
% times into T and V components, reconstruct the mean fields T0, and (vx0,vy0)
% separately and creates two similar struct variables meanField_T and
% meanField_V that are used for the reciprocal reconstruction of the 
% fluctuations T and V. 
% 


path=ttFiltr_r(S,R,xy,tt,Lo,Hi,showFig);
% find reciprocal rays
lLi=length(path);
k=0;
for i=1:lLi-1
    for j=i+1:lLi
        if isequal(path(i).xy0,path(j).xy) && isequal(path(i).xy,path(j).xy0)
            k=k+1;
            in(k,1)=path(i).index;
            t_T(k,1)=(path(i).t+path(j).t)/2;
            t_V(k,1)=(path(i).t-path(j).t)/2;
            Li(k,1)=path(i).length;
            theta(k,:)=path(i).theta;
        end
    end
end
% Temperature reconstruction
d=Li./t_T;
c0est=mean(d);
er=c0est-d;
sig2=er'*er/(k-1);
dc=sqrt(sig2/k);
c02=c0est*c0est;
dc2=2*c0est*dc;
T0est=c02/343/343*293;
dT=dc2/343/343*293;
dt_T=t_T-Li/c0est;
% Velocity reconstruction
d=c02*t_V./Li;
G=-theta;
G_1=pinv(G);
aver=G_1*d;
er=G*aver-d;
sig2=er'*er/(k-2);
std_a=sqrt(sig2*diag(G_1*G_1'));
dt_V=t_V-Li.*(er+d)/c02;
meanField_T.c=c0est;
meanField_T.dc=dc;
meanField_T.T=T0est;
meanField_T.dT=dT;
meanField_V.vx=aver(1);
meanField_V.dvx=std_a(1);
meanField_V.vy=aver(2);
meanField_V.dvy=std_a(2);
meanField_V.c=c0est;
meanField_V.dc=dc;
meanField_V.T=T0est;
meanField_V.dT=dT;
meanField_T.dtt=dt_T;
meanField_T.tt=t_T;
meanField_T.index=in;
meanField_T.data=[-c0est*c0est*dt_T in];
meanField_T.xy=xy;
meanField_T.Li=Li;
meanField_V.dtt=dt_V;
meanField_V.tt=t_V;
meanField_V.index=in;
meanField_V.data=[-c0est*c0est*dt_V in];
meanField_V.xy=xy;
meanField_V.Li=Li;

meanField_T.std_dc=dc;
meanField_T.std_dT=dT;
meanField_T.std_dvx=0;
meanField_T.std_dvy=0;
meanField_V.std_dvx=std_a(1);
meanField_V.std_dvy=std_a(2);
meanField_V.std_dc=0;
meanField_V.std_dT=0;
