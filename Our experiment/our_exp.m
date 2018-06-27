function [field,fields]=our_exp

% load the array
% [file_name,path_name]=uigetfile('*.mat','Locate xy.mat file');
% if file_name==0
%     return
% end
xy=[40.6579 79.7960
    17.5507 1.1259
    77.8724 12.7153
    57.1284 0
    0.8467  67.3008
    0       24.7025
    79.7780 55.0969
    73.0493 79.1313];
xy=xy-40;
% load([path_name,file_name],'xy_Optim_experimental_fit');
%  xy=xy_Optim_experimental_fit_R4_fixed;
%  xy=xy_Optim_experimental_fit;
% xy=xy_Optim;
[file_name,path_name]=uigetfile('*.mat','Pick a file with tt');
if file_name==0
    return
end
load([path_name,file_name],'tt_out');% tt_out is filtered signal tt
% tt_out=zeros(15,40);
% tt_out(1:5,:)=tt(1:5,1:3:end);
% tt_out(6:10,:)=tt(6:10,2:3:end);
% tt_out(11:15,:)=tt(11:15,3:3:end);
[nrays,Ntime]=size(tt_out); %#ok<NODEF>
% noise=0*randn(nrays,1);
% tt_out=tt_out+noise*ones(1,Ntime);
% measured delays of 2.2 ms:
%  tt_out=tt_out-1.2; % increases the reconstructed T to 31 degrees C.
S=3;R=5;
% Ntime=S*Ntime;
tt=zeros(nrays,Ntime);
validInd=1:Ntime;
for i=1:R
   tt(i,:)=interp1(validInd,tt_out(i,:),validInd+1/3,'linear','extrap');
end
tt(R+1:2*R,:)=tt_out(R+1:2*R,:);
for i=2*R+1:3*R
   tt(i,:)=interp1(validInd,tt_out(i,:),validInd-1/3,'linear','extrap');
end
% tt(1:R,1:S:Ntime)=tt_out(1:R,:);
% tt(R+1:2*R,2:S:Ntime)=tt_out(R+1:2*R,:);
% tt(2*R+1:3*R,3:S:Ntime)=tt_out(2*R+1:3*R,:);
tt=tt*1e-3;
% estimation of the mean fields
t0=round(Ntime/2);
Nframes=2;
frames=t0-Nframes:t0+Nframes;
% del_tt=[];
% j0=1;
% for i=1:9
%     j1=j0+4;
%     [meanFields,delays]=spatialMean_serial(S,R,xy,tt(:,j0:j1),320,360,0);
%     del_tt=[del_tt delays(4:end)*1e3];
%     j0=j1;
% end
% del_tt
meanFields=spatialMean_serial(S,R,xy,tt(:,frames),320,360,0,0);
meanFields=[meanFields(Nframes+1) meanFields];
meanFields(1)
[meanFields1,delays]=spatialMean_serial(S,R,xy,tt(:,frames),320,360,1,0);
del_tt=delays(4:end)*1e3 %#ok<NOPRT,NASGU>
meanFields=[meanFields(Nframes+1) meanFields];
% meanFields=spatialMean_serial(S,R,xy,tt,320,360,0);
% meanFields=[meanFields(t0) meanFields(frames)];
% j=0;
% for i=[t0,frames]
%     j=j+1;
%     meanFields(j)=spatialMean(S,R,xy,tt(:,i),320,360,0);
% end

% reconstruction of the fluctuations
Lx=[min(xy(:,1)) max(xy(:,1))];
Ly=[min(xy(:,2)) max(xy(:,2))];
xv=linspace(Lx(1),Lx(2),65);
yv=linspace(Ly(1),Ly(2),65);
% sigma.T=0.3;
% sigma.vx=0.3;
% sigma.vy=0.3;
sigma.T=meanFields(1).std_dT;
sigma.vx=meanFields(1).std_dvx;
sigma.vy=meanFields(1).std_dvy;
sigma.n=5e-6;
sigma.x=1e-2;
lT=20;
lv=20;
%----------------------------------------
U.x=[meanFields.vx]';
U.y=[meanFields.vy]';
% sigma.vx=std(U.x);
% sigma.vy=std(U.y);
U.sig=sqrt((sigma.vx^2+sigma.vy^2)/2);
% sigma.T=0.5*U.sig;
%-------------------------------------------
tau=1.5;
funType='gauss';
estFrame='single';
SI=0;
Cnstr=[];
InvType='Full';
interp=0;
showFig=1;
[field,fields]=t_stochastic_inverse2...
    ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);

