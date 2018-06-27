function [mF,meanFields,meanFields1,meanFields2,d,d1,d2,tt,fluct0,fluct,fluct1,fluct2]=timeDelays
% read the LES
% Lz=[10 12.5];
Lx=[-40 40];Ly=Lx;
S=3;R=5;
load('direct 2D\xy.mat','xy_Optim');
xy=xy_Optim;
sigma.n=5e-6;
sigma.x=0;
wind_dir=3*pi/4;
interp=2;
[tt,VX,VY,~,T,~,T0,vx0,vy0,xv,U,tau,xy]=realFrozenLES(S,R,xy,wind_dir,sigma,1,'vu',1,0,0);
% [T,T0,vx,vx0,vy,vy0,~,~,~,~,xv,yv,~,~,~,~,tau]=LES_3D('u',Lx,Ly,Lz,9:13,0);
pause(1);
% make vy0
% vy0=vy0+2;
% vy=vy+2;
%

lt=size(T,3);
t0=round((lt+1)/2);
T=squeeze(T);
vx=squeeze(VX);
vy=squeeze(VY);
mF.T=T0(t0);
mF.vx=vx0(t0);
mF.vy=vy0(t0);
mF.xy=xy;
mF.index=1:15;
fluct0.dT=T(:,:,t0)-mF.T;
fluct0.vx=vx(:,:,t0)-mF.vx;
fluct0.vy=vy(:,:,t0)-mF.vy;
yv=xv;
menuName='Original fields';
showGraphs(menuName,S,R,mF,xv,yv,fluct0,2,[],[]);

% tt=tt_calc(T,vx,vy,xy,xv,yv,sigma,S,R);
% delays
t1=-1*1e-3;
t2=1*1e-3;
d=rand(S*R,1)*(t2-t1)+t1;
% d([1 10 12])=0;
% d=linspace(t1,t2,15).';
%
%  d=zeros(15,1);
%  d(3)=10*1e-3;
%  d(9)=15*1e-3;
%  d(14)=20*1e-3;
% d=(1:15).'*1e-3;
tt1=tt+d*ones(1,lt);
rel_err=abs(d./tt(:,t0)*100) %#ok<NOPRT>
sort(rel_err)
frames=t0-2:t0+2;
% frames2=t0-4:t0+4;
% tt_tmp=tt1(:,frames);
% tt_tmp(4:end,2)=0;
% c0true=343*sqrt(mean(T0(frames))/293);
% av(1,1)=1/c0true;
% av(2,1)=mean(vx0(frames))/c0true^2;
% av(3,1)=mean(vy0(frames))/c0true^2;
% av(4:18,1)=d;
meanFields=spatialMean_serial(S,R,xy,tt1(:,frames),130,960,0,0);
[meanFields1,d1]=spatialMean_serial(S,R,xy,tt1(:,frames),130,960,1,0);
% [meanFields2,d2]=spatialMean_serial_t(S,R,xy,tt1(:,frames),130,960,0);
[meanFields2,d2]=spatialMean_serial_maxD(S,R,xy,tt1(:,frames),130,960,1,sum(d~=0),0);
% [meanFields2,d2]=spatialMean_serial_seq(S,R,xy,tt1(:,frames),130,960,1,0);
table1=[mean(T0(frames))-273.15 mean(vx0(frames)) mean(vy0(frames)); ...
    meanFields(1).T-273.15 meanFields(1).vx meanFields(1).vy; ...
    meanFields(1).std_dT meanFields(1).std_dvx meanFields(1).std_dvy;...
    meanFields1(1).T-273.15 meanFields1(1).vx meanFields1(1).vy;...
    meanFields1(1).std_dT meanFields1(1).std_dvx meanFields1(1).std_dvy;...
    meanFields2(1).T-273.15 meanFields2(1).vx meanFields2(1).vy;...
    meanFields2(1).std_dT meanFields2(1).std_dvx meanFields2(1).std_dvy] %#ok<NASGU,NOPRT>
figure;plot(tt(:,t0)*1000,'*','linewidth',2,'MarkerSize',10);
xlabel('Ray index','fontsize',12,'fontweight','bold')
ylabel('Travel time (ms)','fontsize',12,'fontweight','bold')
figure;plot(d*1000,'b*','linewidth',2,'MarkerSize',10);
hold on
plot(d1*1000,'o','linewidth',2,'MarkerSize',10,'Color',[0 0.5 0]);
plot(d2*1000,'rs','linewidth',2,'MarkerSize',10);

xlabel('Ray index','fontsize',12,'fontweight','bold')
ylabel('Systematic error (ms)','fontsize',12,'fontweight','bold')
legend('Actual','Estimated','Estimated d_t',0);
% fluctuations
% Lx=[xv(1) xv(end)];
% Ly=[yv(1) yv(end)];
sigma.T=0.04;
sigma.vx=0.03;
sigma.vy=0.03;

% U.x=4;U.y=0;U.sig=sigma.vx*sqrt(2);
% tau=1;
lT=16;lv=16;

funType='Gauss';
estFrame='single';
Cnstr=[];
showFig=1;
SI=0;
invType='Full';
t0index=(1+length(frames))/2;
fluct=t_stochastic_inverse2('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    [meanFields(t0index) meanFields],U,tau,t0,frames,funType,estFrame,SI,Cnstr,invType,interp,showFig);
fluct1=t_stochastic_inverse2('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    [meanFields1(t0index) meanFields1],U,tau,t0,frames,funType,estFrame,SI,Cnstr,invType,interp,showFig);
fluct2=t_stochastic_inverse2('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    [meanFields2(t0index) meanFields2],U,tau,t0,frames,funType,estFrame,SI,Cnstr,invType,interp,showFig);

