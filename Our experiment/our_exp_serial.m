function [field,field1,field2,field3,meanFields,meanFields1,meanFields2,meanFields3]=our_exp_serial

% speakers
xyzS=[4.0726200e+01   7.9458700e+01   2.4450000e-01
   7.2675700e+01   7.8965000e+01   2.4889000e-01
   7.9565200e+01   5.5355500e+01   2.3424000e-01
   7.7575400e+01   1.2851300e+01   2.4780000e-01
   5.7387400e+01   2.1480000e-01   2.3877000e-01
   1.7644800e+01   1.4590000e+00   2.7250000e-01
   2.4800000e-01   2.4487600e+01   2.2962000e-01
   1.0125000e+00   6.6921600e+01   2.5722000e-01];
% mics
xyzM=[4.1245300e+01   7.9486400e+01   2.4450000e-01
   7.3411600e+01   7.8837300e+01   2.4249000e-01
   7.9504300e+01   5.4732400e+01   2.3424000e-01
   7.7491200e+01   1.2140700e+01   2.4780000e-01
   5.6847200e+01   1.3970000e-01   2.2927000e-01
   1.7028500e+01   1.5200000e+00   2.7250000e-01
   8.8900000e-02   2.4959400e+01   2.3282000e-01
   1.1121000e+00   6.7602500e+01   2.3492000e-01];
% final 2D aray
xy=[xyzS(:,1:2);xyzM(:,1:2)];
xy=xy-40;
S=8;
R=8;
% get a file with tt
[file_name,path_name]=uigetfile('*.mat','Pick a file with tt');
if file_name==0
    return
end
%load([path_name,file_name],'tt_out'); % tt_out is filtered signal tt
load([path_name,file_name],'tt'); % tt_out is filtered signal tt
%% compencation for hardware time delays
% delay times in hardware (ms): [i,j] stands for the i-th Speaker and j-th
% Microphone
tDel=[1.5659994e+00   1.5695211e+00   1.6635981e+00   1.7860622e+00   1.7934960e+00   1.7604294e+00   1.6716353e+00   1.5691098e+00
   1.6447130e+00   1.6716521e+00   1.7625810e+00   1.8171108e+00   1.8532592e+00   1.8545286e+00   1.7919826e+00   1.6791170e+00
   1.7555482e+00   1.5956268e+00   1.5790157e+00   1.5766389e+00   1.6013185e+00   1.7092368e+00   1.7965629e+00   1.8188520e+00
   1.6989301e+00   1.5965757e+00   1.5707215e+00   1.7436454e+00   1.5867701e+00   1.6729609e+00   1.7588060e+00   1.7727558e+00
   1.7847536e+00   1.7651335e+00   1.7140544e+00   1.5693702e+00   1.5429242e+00   1.5604774e+00   1.6207589e+00   1.7350560e+00
   1.7893890e+00   1.7749554e+00   1.7128875e+00   1.5901768e+00   1.5596270e+00   1.5834358e+00   1.6388312e+00   1.7374633e+00
   1.7752748e+00   1.8425482e+00   1.8468503e+00   1.7955835e+00   1.7284468e+00   1.6263876e+00   1.6102992e+00   1.6492032e+00
   1.7324933e+00   1.7709918e+00   1.8119330e+00   1.8029193e+00   1.7351612e+00   1.6327652e+00   1.5928638e+00   1.6012841e+00];
delT=tDel.';
delT=delT(:);
[nrays,Ntime]=size(tt); %#ok<NODEF>
tt=tt-repmat(delT,1,Ntime);
tt=tt*1e-3; % convert in s
%%
%---------------------------------------
% estimation of the mean fields
meanField=[];
meanField1=[];
meanFields3=[];
N_frames=2;
% mean_frames=3;
% mean_frames=2*N_frames+1; % used for the mean field reconstruction
mean_frames=Ntime;
if Ntime/mean_frames~=fix(Ntime/mean_frames)
    error('Ntime=%g, mean_frames=%g, Ntime/mean_frames must be integer',Ntime,mean_frames);
end
for i=1:Ntime/mean_frames
    tmp=spatialMean_serial(S,R,xy,tt(:,(i-1)*mean_frames+1:i*mean_frames),320,360,0,0);
    meanField=[meanField tmp]; %#ok<AGROW>
    [tmp,delays]=spatialMean_serial(S,R,xy,tt(:,(i-1)*mean_frames+1:i*mean_frames),320,360,1,0);
%     [tmp,delays]=spatialMean_serial_maxD(S,R,xy,tt(:,(i-1)*mean_frames+1:i*mean_frames),320,360,1,15,0);
    meanField1=[meanField1 tmp]; %#ok<AGROW>
%     [tmp,delays3]=spatialMean_serial_t(S,R,xy,tt(:,(i-1)*mean_frames+1:i*mean_frames),320,360,0);
%     meanField3=[meanFields3 tmp]; %#ok<AGROW>
    figure(i);
    plot(delays*1000,'o','linewidth',2,'color',[0 0.5 0],'MarkerSize',10)
    xlabel('Path index','fontsize',12,'fontweight','bold')
    ylabel('Time delay (ms)','fontsize',12,'fontweight','bold')
    
end
tmp=input('i0  and Ndelays: ');
i0=tmp(1);
Ndel=tmp(2);
t0=i0*mean_frames-round((mean_frames-1)/2);
% t0=round(Ntime/2);
frames=t0-N_frames:t0+N_frames;
meanFields=[meanField(t0) meanField(frames)];
meanFields1=[meanField1(t0) meanField1(frames)];
disp('%----------------------------------------')
disp(meanFields(1));
disp('%----------------------------------------')
disp(meanFields1(1));
disp('%----------------------------------------')
if Ndel~=0 % delays numbers are indicated
    delaysFlg=1;
else
    delaysFlg=0; % no delays are indicated
end
[tmp,d2]=spatialMean_serial_maxD(S,R,xy,tt(:,(i0-1)*mean_frames+1:i0*mean_frames),320,360,delaysFlg,Ndel,0);
meanFields2=[tmp(mean_frames-N_frames) tmp];
disp(meanFields2(1));
disp('%----------------------------------------')
figure(i0)
hold on
plot(d2*1000,'rs','linewidth',2,'MarkerSize',10);

% reconstruction of the fluctuations
Lx=[min(xy(:,1)) max(xy(:,1))];
Ly=[min(xy(:,2)) max(xy(:,2))];
xv=linspace(Lx(1),Lx(2),65);
yv=linspace(Ly(1),Ly(2),65);
% sigma.T=0.04;
% sigma.vx=0.04;
% sigma.vy=0.04;
% sigma.T=0.2;
% sigma.vx=0.3;
% sigma.vy=0.3;
sigma.T=meanFields(1).std_dT;
sigma.vx=meanFields(1).std_dvx;
sigma.vy=meanFields(1).std_dvy;

sigma.n=1e-5;
sigma.x=1e-3;
lT=15;
lv=15;
%----------------------------------------
U.x=[meanFields.vx]';
U.y=[meanFields.vy]';
% sigma.vx=std(U.x);
% sigma.vy=std(U.y);
U.sig=sqrt((sigma.vx^2+sigma.vy^2)/2);
% sigma.T=0.5*U.sig;
U1.x=[meanFields1.vx]';
U1.y=[meanFields1.vy]';
U1.sig=U.sig;
sigma1=sigma;
sigma1.T=meanFields1(1).std_dT;
sigma1.vx=meanFields1(1).std_dvx;
sigma1.vy=meanFields1(1).std_dvy;

U2.x=[meanFields2.vx]';
U2.y=[meanFields2.vy]';
U2.sig=U.sig;
sigma2=sigma;
sigma2.T=meanFields2(1).std_dT;
sigma2.vx=meanFields2(1).std_dvx;
sigma2.vy=meanFields2(1).std_dvy;

%-------------------------------------------
tau=0.5;
funType='gauss';
estFrame='single';
SI=0;
Cnstr=[];
InvType='Full';
interp=0;
showFig=1;

% [field,~,~,R_dd0, R_md]=t_stochastic_inverse2...
%     ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
%     meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);

field1=t_stochastic_inverse2...
    ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma1,lT,lv,...
    meanFields1,U1,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);

field2=t_stochastic_inverse2...
    ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma2,lT,lv,...
    meanFields2,U2,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);

