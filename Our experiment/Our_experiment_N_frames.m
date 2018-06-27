function Our_experiment_N_frames(Ntime,lT,lv)
% Our_experiment_N_frames solves the 2D tomography problem for the STINHO experiment.
% The fields are reconstructed for a time interval [11:44-12:16].
% Ntime consecutive frames are used to reconstruct the middle frame.
% Syntax:
% Our_experiment_N_frames(Ntime,lT,lv); or Our_experiment_N_frames;
% In the latter case, you will be promted to load the previously saved file (see output description).  
% Inputs:
% Ntime - the number of consecutive frames to be used for the
%   reconstruction of a single frame, an odd integer;
% lT, lv - correlation lengths of T and v fluctuations, meters;
% Output:
% a file 'German_frames_Ntime_lT_lv.mat' (with corresponding numerical values of Ntime, lT, and lv)
% will be created. This file will contain the following variables:
% mF - reconstructed mean fields, a struct vector with fields:
%   mF.c = reconstructed c0, a scalar;
%   mF.T = reconstructed T0, a scalar; 
%   mF.vx = reconstructed vx0, a scalar;
%   mF.vy = reconstructed vy0, a scalar;
%   mF.std_dc = estimation of the errors in the reconstruction of c0, a scalar;
%   mF.std_dT = estimation of the errors in the reconstruction of T0, a scalar;
%   mF.std_dvx = estimation of the errors in the reconstruction of vx0, a scalar;
%   mF.std_dvy = estimation of the errors in the reconstruction of vy0, a scalar;
%   mF.dtt = travel times in seconds due to fluctuations only;
%   mF.tt = full travel times used as the input;
%   mF.index = indeces of valid travel times, a vector of integers;
%   mF.data = input data for the reconstruction of the fluctuations and their indeces
%       a matrix [data index], where data = -c0^2*dtt;
%   mF.xy = xy input array;
%   mF.Li = the lengths of the valid travel paths, a vector of the length less or equal
%       to S*R (can be less if some of the travel times are filtered out);
% fluct - a struct variable with fields:
%   fluct.dT = T fluctuations, a matrix (x,y,t);
%   fluct.vx = vx fluctuations, a matrix (x,y,t);
%   fluct.vy = vy fluctuations, a matrix (x,y,t);
%   fluct.dC = dC fluctuations, a matrix (x,y,t);
%   fluct.dTExpNMSE = dT expected NMSE, a matrix (x,y,t);
%   fluct.dTExpAverNMSE = dT spatial mean of expected NMSE, a vector (1,t);
%   fluct.dTExpSTD = dT expected STD, a matrix (x,y,t);
%   fluct.dTExpAverSTD = dT spatial mean of expected STD, a vector (1,t);
%   fluct.vxExpNMSE = vx expected NMSE, a matrix (x,y,t);
%   fluct.vxExpAverNMSE = vx spatial mean of expected NMSE, a vector (1,t);
%   fluct.vxExpSTD = vx expected STD, a matrix (x,y,t);
%   fluct.vxExpAverSTD = vx spatial mean of expected STD, a vector (1,t);
%   fluct.vyExpNMSE = vy expected NMSE, a matrix (x,y,t);
%   fluct.vyExpAverNMSE = vy spatial mean of expected NMSE, a vector (1,t);
%   fluct.vyExpSTD = vy expected STD, a matrix (x,y,t);
%   fluct.vyExpAverSTD = vy spatial mean of expected STD, a vector (1,t);
% fluct2 - the same as fluct, but only for the middle frame of the whole time interval.
% Example:
% Our_experiment_N_frames(11,20,20);

% field=[];
% fields=[];
if nargin==3
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
Lx=[min(xy(:,1)) max(xy(:,1))];
Ly=[min(xy(:,2)) max(xy(:,2))];
xv=linspace(Lx(1),Lx(2),65);
yv=linspace(Ly(1),Ly(2),65);
    Lo=320;
    Hi=360;
%     sigma.T=0.04;
%     sigma.vx=0.03;
%     sigma.vy=0.03;
%     sigma.T=0.1;
%     sigma.vx=0.1;
%     sigma.vy=0.1;

%     sigma.n=1e-5;
%     sigma.x=1e-3;
%     U.sig=sigma.vx*sqrt(2);
    tau=0.5;
    funType='gauss';
    estFrame='single';
    Cnstr=[];
    interp=0;
    showFig=0;
    SI=0;
    invType='Full';
    %% read data files
    t=[];
    D=dir('tt*.mat');
    fileInd=input(['There are ',num2str(length(D)),' files with the data in the current dir. Enter indices which ones to use: ']);
    if ~isempty(fileInd)
        D=D(fileInd);
    end
    for i=1:length(D)
        try
            load(D(i).name,'tt_out');
            % compatibility with old file format
            if isequal('tt_out',D(i).name(1:6))
                t=tt_out;
                break
            else
                t=[t tt_out];
            end
        catch
            load(D(i).name,'tt');
            t=[t tt];
        end
    end
    [~,l_tt]=size(t);
    
% Determine the time of the experiment
Year=D(1).name(4:7);
Month=D(1).name(8:9);
Day=D(1).name(10:11);
StartTime=[D(1).name(12:13),':',D(1).name(14:15)];
EndTime=[D(end).name(12:13),':',D(end).name(14:15)];

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

%%
t=t-repmat(delT,1,l_tt); % correction for time delays
%
mkdir('tt');
for i=1:S*R
    figure(1);
    plot((0:l_tt-1)*tau,t(i,:),'o');
    xlabel('Time (s)','fontweight','bold','fontsize',12);
    ylabel('Travel time (ms)','fontweight','bold','fontsize',12);
    title(sprintf('%s-%s-%s, %s - %s, Ray %d',Month,Day,Year,StartTime,EndTime,i),'fontweight','bold','fontsize',12);
    saveas(gcf,sprintf('tt\\%s-%s-%s Ray %d',Month,Day,Year,i),'fig');
end
pause(0.5);
t=t/1000; % in sec
%% calibrate the array: determine effective xy
%  xy=arrayCalibration(t,xy);
%     Lx=[min(xy(:,1)) max(xy(:,1))];
%     Ly=[min(xy(:,2)) max(xy(:,2))];
%     xv=linspace(Lx(1),Lx(2),101);
%     yv=linspace(Ly(1),Ly(2),101);

% determine U and sigmas
tmp=spatialMean_serial(S,R,xy,t,Lo,Hi,0,0);
U.x=mean([tmp.vx]);
U.y=mean([tmp.vy]);
% sigma.T=0.2;
% sigma.vx=0.2;
% sigma.vy=0.2;
sigma.T=tmp(1).std_dT;
sigma.vx=tmp(1).std_dvx;
sigma.vy=tmp(1).std_dvy;

sigma.n=1e-5;
sigma.x=1e-3;
U.sig=sqrt(tmp(1).std_dvx^2+tmp(1).std_dvy^2);
%% Run the reconstruction
hw = waitbar(0,'Please wait...');
     Nframes=(Ntime-1)/2;
     for ij=1:l_tt-Ntime+1
        t0=Nframes+ij;
        frames=t0-Nframes:t0+Nframes;
        i=0;
        for j=[t0,frames]
            i=i+1;
            meanFields(i)=spatialMean_serial(S,R,xy,t(:,j),Lo,Hi,0,0);
        end
%         meanFields=spatialMean_serial(S,R,xy,t(:,frames),Lo,Hi,0,0);
%         meanFields=[meanFields(Nframes+1) meanFields];
%      Syntax:  (ddType,R_dd0,mdType,R_md,S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
%              meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);
        if ij==1
            [field,~,~,R_dd0, R_md]=t_stochastic_inverse2('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
           meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,invType,interp,showFig);
%              load(mat_file,'R_dd0','R_md');
        else
            field=t_stochastic_inverse2('o',R_dd0,'o',R_md,S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
           meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,invType,interp,showFig);
        end
        mF(ij)=meanFields(1);
        fluct.dT(:,:,ij)=field.dT;
        fluct.vx(:,:,ij)=field.vx;
        fluct.vy(:,:,ij)=field.vy;
        fluct.dC(:,:,ij)=field.dC;
        fluct.dTExpNMSE(:,:,ij)=field.dTExpNMSE;
        fluct.dTExpAverNMSE(ij)=field.dTExpAverNMSE;
        fluct.dTExpSTD(:,:,ij)=field.dTExpSTD;
        fluct.dTExpAverSTD(ij)=field.dTExpAverSTD;
        fluct.vxExpNMSE(:,:,ij)=field.vxExpNMSE;
        fluct.vxExpAverNMSE(ij)=field.vxExpAverNMSE;
        fluct.vxExpSTD(:,:,ij)=field.vxExpSTD;
        fluct.vxExpAverSTD(ij)=field.vxExpAverSTD;
        fluct.vyExpNMSE(:,:,ij)=field.vyExpNMSE;
        fluct.vyExpAverNMSE(ij)=field.vyExpAverNMSE;
        fluct.vyExpSTD(:,:,ij)=field.vyExpSTD;
        fluct.vyExpAverSTD(ij)=field.vyExpAverSTD;
        waitbar(ij/(l_tt-Ntime+1),hw);
    end %ij
    close(hw);
        mF1=[];fluct1=[];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10-min averaged fields
%     mF1=mF(1);
%     N=length(mF);
%     mF1.T=mean([mF.T]);
%     mF1.vx=mean([mF.vx]);
%     mF1.vy=mean([mF.vy]);
%     tmp=mean([mF.std_dT].^2)/N;
%     mF1.std_dT=sqrt(tmp);
%     tmp=mean([mF.std_dvx].^2)/N;
%     mF1.std_dvx=sqrt(tmp);
%     tmp=mean([mF.std_dvy].^2)/N;
%     mF1.std_dvy=sqrt(tmp);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fluct1.dT=mean(fluct.dT,3);
%     fluct1.vx=mean(fluct.vx,3);
%     fluct1.vy=mean(fluct.vy,3);
%     fluct1.dC=mean(fluct.dC,3);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tmp=mean(fluct.dTExpSTD.^2,3)/N;
%     sigmas=(fluct.dTExpSTD.^2)./fluct.dTExpNMSE;
%     m_sig=mean(sigmas,3);
%     fluct1.dTExpNMSE=tmp./m_sig;
%     fluct1.dTExpAverNMSE=mean(fluct1.dTExpNMSE(:));
%     fluct1.dTExpSTD=sqrt(tmp);
%     fluct1.dTExpAverSTD=mean(fluct1.dTExpSTD(:));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tmp=mean(fluct.vxExpSTD.^2,3)/N;
%     sigmas=(fluct.vxExpSTD.^2)./fluct.vxExpNMSE;
%     m_sig=mean(sigmas,3);
%     fluct1.vxExpNMSE=tmp./m_sig;
%     fluct1.vxExpAverNMSE=mean(fluct1.vxExpNMSE(:));
%     fluct1.vxExpSTD=sqrt(tmp);
%     fluct1.vxExpAverSTD=mean(fluct1.vxExpSTD(:));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tmp=mean(fluct.vyExpSTD.^2,3)/N;
%     sigmas=(fluct.vyExpSTD.^2)./fluct.vyExpNMSE;
%     m_sig=mean(sigmas,3);
%     fluct1.vyExpNMSE=tmp./m_sig;
%     fluct1.vyExpAverNMSE=mean(fluct1.vyExpNMSE(:));
%     fluct1.vyExpSTD=sqrt(tmp);
%     fluct1.vyExpAverSTD=mean(fluct1.vyExpSTD(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%% Fields at the middle point of time interval
    fr=round((l_tt-Ntime+2)/2);
    fluct2.dT=fluct.dT(:,:,fr);
    fluct2.vx=fluct.vx(:,:,fr);
    fluct2.vy=fluct.vy(:,:,fr);
    fluct2.dC=fluct.dC(:,:,fr);
    fluct2.dTExpNMSE=fluct.dTExpNMSE(:,:,fr);
    fluct2.dTExpAverNMSE=fluct.dTExpAverNMSE(fr);
    fluct2.dTExpSTD=fluct.dTExpSTD(:,:,fr);
    fluct2.dTExpAverSTD=mean(fluct2.dTExpSTD(:));
    fluct2.vxExpNMSE=fluct.vxExpNMSE(:,:,fr);
    fluct2.vxExpAverNMSE=mean(fluct2.vxExpNMSE(:));
    fluct2.vxExpSTD=fluct.vxExpSTD(:,:,fr);
    fluct2.vxExpAverSTD=mean(fluct2.vxExpSTD(:));
    fluct2.vyExpNMSE=fluct.vyExpNMSE(:,:,fr);
    fluct2.vyExpAverNMSE=mean(fluct2.vyExpNMSE(:));
    fluct2.vyExpSTD=fluct.vyExpSTD(:,:,fr);
    fluct2.vyExpAverSTD=mean(fluct2.vyExpSTD(:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load(mat_file,'SingValG','SingValC');
    SingValG=[];SingValC=[];
    fileName=sprintf('OurExp_frames_%d_%g_%g.mat',Ntime,lT,lv);
    menuName1=['Fields at ',StartTime,'-',EndTime];
    SingleTime=getFrameTime(fr,StartTime,tau);
    menuName=['Fields at ',SingleTime];
    for i=1:l_tt-Ntime+1
        dt(:,i)=mF(i).dtt;
    end
    mkdir('dtt');
    for i=1:size(dt,1)
        figure(1);
        plot((0:l_tt-Ntime)*tau,dt(i,:)*1000,'-o');
        xlabel('Time (s)','fontweight','bold','fontsize',12);
        ylabel('Travel time differences (ms)','fontweight','bold','fontsize',12);
        title(sprintf('%s-%s-%s, %s - %s, Ray %d',Month,Day,Year,StartTime,EndTime,i),'fontweight','bold','fontsize',12);
        saveas(gcf,sprintf('dtt\\%s-%s-%s Ray %d',Month,Day,Year,i),'fig');
    end
    mF(1).t=(0:l_tt-Ntime)*tau; % time axis for mean field plots
    
    save (fileName,'fluct1','mF1','fluct','fluct2','fr','mF','S','R','xv','yv','SingValG','SingValC','menuName','menuName1','t','dt');
    % Syntax:
    % showGraphs(menuName,S,R,meanFields,xv,yv,fluct,interp,SingValG,SingValC)
  
    showGraphs(menuName1,S,R,mF,xv,yv,fluct,interp,SingValG,SingValC);
     
    showGraphs(menuName,S,R,mF(fr),xv,yv,fluct2,interp,SingValG,SingValC);

elseif nargin==0
    [FileName,PathName] = uigetfile('*.mat','Select a MAT-file');
    if FileName~=0
        load([PathName,FileName]);
    end
    interp=input('Interpolation factor: ');
%     menuName='Fields at 22:15-22:21';
    showGraphs(menuName1,S,R,mF,xv,yv,fluct,interp,SingValG,SingValC);
%     menuName='Average fields at 5:30 AM';
%     showGraphs(menuName,S,R,mF1,xv,yv,fluct1,interp,SingValG,SingValC);
%     menuName='Fields at 22:19';
    showGraphs(menuName,S,R,mF(fr),xv,yv,fluct2,interp,SingValG,SingValC);
    
end

function SingleTime=getFrameTime(fr,StartTime,tau)
elapsedTime=fr*tau; % in sec
% convert StartTime in sec
h0=str2num(StartTime(1:2))*3600; %#ok<*ST2NM>
m0=str2num(StartTime(4:5))*60;
t0=h0+m0;% StartTime in sec
% add the elapsed time
t=t0+elapsedTime;
% convert it back
h=fix(t/3600); % assumed to be not over 24 hours (same day)
m=fix((t-h*3600)/60);
s=round(t-h*3600-m*60);
% insert 0 in front, if needed
h=num2str(h,'%d');
if length(h)==1
    h=['0',h];
end
m=num2str(m,'%d');
if length(m)==1
    m=['0',m];
end
s=num2str(s,'%d');
if length(s)==1
    s=['0',s];
end
% parse it to hours, min, and sec
SingleTime=[h,':',m,':',s];



