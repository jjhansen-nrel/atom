function meanField=spatialMean(S,R,xy,tt,Lo,Hi,showFig)
% spatialMean filters the initial travel times, reconstructs the mean
% fields T0, vx0, vy0 in the 2D inverse problem, estimates the errors
% of this recostruction and forms a struct
% variable meanField which is used later in TDSI (t_stochastic_inverse2).
% Syntax:
% meanField=spatialMean(S,R,xy,tt,Lo,Hi,showFig);
% Inputs:
% S - number of sound sources, integer;
% R - number of receivers, integer;
% xy - the x and y coordinated of sources and receivers; a
%   matrix [S+R, 2]; the first S rows contain the coordinates of sources,
%   the last R rows contain the coordinates of receivers;
% tt - travel times (sec), a vector of length S*R;
% Lo - a lower limit for the filtration; travel times which corresponds to
%   the effective sound speed lower then Lo will be ommited;
% Hi - an upper limit for the filtration; travel times which corresponds to
%   the effective sound speed higher then Hi will be ommited;
% showFig - a flag to display a graph; if showFig=1, the graph showing Lo,
%   Hi and effective sound speed will be displayed, otherwise, it will not be
%   displayed.
% Outputs:
% meanField - a struct variable with fields:
%   meanField.c = reconstructed c0, a scalar;
%   meanField.T = reconstructed T0, a scalar; 
%   meanField.vx = reconstructed vx0, a scalar;
%   meanField.vy = reconstructed vy0, a scalar;
%   meanField.std_dc = estimation of the errors in the reconstruction of c0, a scalar;
%   meanField.std_dT = estimation of the errors in the reconstruction of T0, a scalar;
%   meanField.std_dvx = estimation of the errors in the reconstruction of vx0, a scalar;
%   meanField.std_dvy = estimation of the errors in the reconstruction of vy0, a scalar;
%   meanField.dtt = filtered travel times in seconds due to fluctuations only, 
%       a vector of size [Number of valid travel times, 1];
%   meanField.tt = filtered full travel times, a vector of the same size as dtt;
%   meanField.index = indeces of valid travel times, a vector of integers of the same size as dtt;
%   meanField.data = input data for the reconstruction of the fluctuations and their indeces
%   a matrix [data index], where data = -c0^2*dtt;
%   meanField.xy = xy array of transducers;
%   meanField.Li = the lengths of the valid travel paths, a vector of the length less or equal
%   to S*R (can be less if some of the travel times are filtered out);

[Li,tt,in,theta]=ttFiltr(S,R,xy,tt,Lo,Hi,showFig);
lLi=length(Li);
G=[ones(lLi,1) -theta];
d=tt./Li;
G_1=pinv(G);
aver=G_1*d;
er=G*aver-d;
sig2=er'*er/(lLi-3);
std_a=sqrt(sig2*diag(G_1*G_1'));
c0est=1/aver(1);
c02=c0est*c0est;
dc=c02*std_a(1);
dc2=2*c0est*dc;
gamma=1.4;
Ra=287.058;
gR=gamma*Ra;
T0est=c02/gR;
dT=dc2/gR;
vx0est=aver(2)*c02;
dvx=c02*std_a(2);%+aver(2)*dc2;
vy0est=aver(3)*c02;
dvy=c02*std_a(3);%+aver(3)*dc2;
tt0=Li.*(er+d);% travel times due to the mean fields only
dtt=tt-tt0;
meanField.c=c0est;
%meanField.dc=dc;
meanField.T=T0est;
%meanField.dT=dT;
meanField.vx=vx0est;
%meanField.dvx=dvx;
meanField.vy=vy0est;
%meanField.dvy=dvy;
%%%%%%%%%%%%%%%%%%%%%%%% estimation of sigma due to turbulence
% sig_tt=sig.tt;
% sig_x=sig.xy;
% var_tt_x=(sig_tt./Li).^2+(d./Li.*dLi).^2;
% if max(var_tt_x)>sig2
%     error('The variance of the fluctuations cannot be negative! Make errors of tt or xy smaller.');
% end
% 
% cov_d=diag(sig2*ones(lLi,1)-var_tt_x);
% std_a=sqrt(diag(G_1*cov_d*G_1'));
% dc=c02*std_a(1)
% dc2=2*c0est*dc;
% dT=dc2/343/343*293
% dvx=c02*std_a(2)+aver(2)*dc2
% dvy=c02*std_a(3)+aver(3)*dc2
% 
% % cov_d=sig2*ones(lLi,1)-var_tt_x;
% % for i=1:lLi
% %     Y(i,1)=quadl0('gaussCorr2D_T',0,Li(i),1e-8,[],path(in(i)),c0est,T0est,lT);
% %     Y(i,2)=quadl0('gaussCorr2D_V',0,Li(i),1e-8,[],path(in(i)),lv);
% %     Y(i,:)=Y(i,:)/Li(i)^2;
% % end
% % sTV2=pinv(Y)*cov_d;
% % dT=sqrt(sTV2(1));
% % dvx=sqrt(sTV2(2));
% % dvy=dvx;
% % dc=343*343/(293*2*c0est)*dT;
% 
% meanField.sig_dc=dc;
% meanField.sig_dT=dT;
% meanField.sig_dvx=dvx;
% meanField.sig_dvy=dvy;
% %%%%%%%%%%%%%%%%%%%%%%%% estimation of sigma due to travel time and positions errors
% cov_d=diag(var_tt_x);
% std_a=sqrt(diag(G_1*cov_d*G_1'));
% dc=c02*std_a(1);
% dc2=2*c0est*dc;
% dT=dc2/343/343*293;
% dvx=c02*std_a(2)+aver(2)*dc2;
% dvy=c02*std_a(3)+aver(3)*dc2;
meanField.std_dc=dc;
meanField.std_dT=dT;
meanField.std_dvx=dvx;
meanField.std_dvy=dvy;
%%%%%%%%%%%%%%%%%%%%%%%%% misc
% meanField.std_tt=sig_tt;
% meanField.std_x=sig_x;
meanField.dtt=dtt;
meanField.tt=tt;
meanField.index=in;
meanField.data=[-c0est*c0est*dtt in];
meanField.xy=xy;
meanField.Li=Li;
meanField.tt0=tt0;