function [meanField,field,mat_file]=LESinverse(S,R,xy,lT,lv,funType,FigInterp,CnstrFlag,sigma_t,sigma_x)
% LESinverse reads LES data, solves the forward problem for them and solves the inverse problem.
% Syntax:
% [meanField,field,mat_file]=LESinverse(S,R,xy,lT,lv,funType,FigInterp,CnstrFlag,sigma_t,sigma_x);
% Inputs:
% S - number of sound sources, integer;
% R - number of receivers, integer;
% xy - x and y coordinates of sorces and receivers; xy is a matrix of size
%   [S+R,2]; each row of xy contains x and y coordinates of one transducer;
%   the first S rows should contain the coordinates of S sources, the last R
%   rows contain the coordinates of the receivers; if xy=[] (empty array),
%   the optimal location will be loaded from file 'direct 2D\xy.mat';
% lT, lv - correlation lengths of T and v fluctuations, meters;
% funType - a string which specifies what covariance functions should be
%   used in calculations; the optiona are:
%   'gauss' - Gaussian covariance functions;
%   'exp' - exponent;
%   'exp23' - exp[-(R/L)^(2/3)];
%   'vonk' - von Karman;
%   The calculation is the fastest with 'gauss' type;
% FigInterp - interpolation of the displayed fields for illustrative purposes;
%   FigInterp = 0 - no interpolation; FigInterp = 1 - one extra point will be added
%   in the middle of two adjacent points in space and time, ant etc (step=1/2^FigInterp).
% CnstrFlag - if CnstrFlag = 1 additional linear constraints will be imposed
%   on the fields of fluctuations: Laplacian of T and div of V; the values of
%   these operators will be calculated bhy the finite difference method using
%   the original fields; if CnstrFlag = 0, no constraints are imposed.
% sigma_t - the standard deviation (in seconds) of the errors in the travel times;
% sigma_x - the standard deviation (in meters) of the errors in the transducer locations;
% Outputs:
% meanField - reconstruction of mean fields at specified time t0. meanField has the following fields:
%   meanField.c = reconstructed c0, a scalar;
%   meanField.T = reconstructed T0, a scalar; 
%   meanField.vx = reconstructed vx0, a scalar;
%   meanField.vy = reconstructed vy0, a scalar;
%   meanField.std_dc = estimation of the errors in the reconstruction of c0, a scalar;
%   meanField.std_dT = estimation of the errors in the reconstruction of T0, a scalar;
%   meanField.std_dvx = estimation of the errors in the reconstruction of vx0, a scalar;
%   meanField.std_dvy = estimation of the errors in the reconstruction of vy0, a scalar;
%   meanField.dtt = filtered travel times in seconds due to fluctuations only;
%   meanField.tt = filtered full travel times used as the input for mean field reconstruction;
%   meanField.index = indeces of valid travel times, a vector of integers;
%   meanField.data = input data for the reconstruction of the fluctuations and their indeces
%   a matrix [data index], where data = -c0^2*dtt;
%   meanField.xy = xy input array;
%   meanField.Li = the lengths of the valid travel paths, a vector of the length less or equal
%   to S*R (can be less if some of the travel times are filtered out);
% field - the reconstructed fluctuations at time t0, a struct variable with fields:
%   field.dT = T fluctuations, a matrix (x,y);
%   field.vx = vx fluctuations, a matrix (x,y);
%   field.vy = vy fluctuations, a matrix (x,y);
%   field.dC = dC fluctuations, a matrix (x,y);
%   field.dTExpNMSE = dT expected NMSE, a matrix (x,y);
%   field.dTExpAverNMSE = dT spatial mean of expected NMSE, a scalar;
%   field.dTExpSTD = dT expected STD, a matrix (x,y);
%   field.dTExpAverSTD = dT spatial mean of expected STD, a scalar;
%   field.vxExpNMSE = vx expected NMSE, a matrix (x,y);
%   field.vxExpAverNMSE = vx spatial mean of expected NMSE, a scalar;
%   field.vxExpSTD = vx expected STD, a matrix (x,y);
%   field.vxExpAverSTD = vx spatial mean of expected STD, a scalar;
%   field.vyExpNMSE = vy expected NMSE, a matrix (x,y);
%   field.vyExpAverNMSE = vy spatial mean of expected NMSE, a scalar;
%   field.vyExpSTD = vy expected STD, a matrix (x,y);
%   field.vyExpAverSTD = vy spatial mean of expected STD, a scalar;
% mat_file - a name of '.mat' file with the results; this file is
%   atomatically created in a subdirectory 'TDSI_Res' in the current dirctory
%   (if there is no such a subdirectory, it will be created); the name for the
%   first run of t_stochastic_inverse2 creates file '1.mat', the second -
%   '2.mat', and etc. After the 20th run, the file names will start wih 1
%   again, so that in the direcoty 'TDSI_Res' will be stored the results of
%   20 recent runs of t_stochastic_inverse2; in this mat_file all vaiables
%   are stored, including constrained solution 'field_Cnstr' and 'fields_Cnstr',
%   and ordinary stochastic inversion 'field_SI' and 'fields_SI'; these variables
%   are struct variables with the same fields as 'field'; for each 
%   '.mat' file, a '.txt' file is also created with the same name containing a 
%   short description of the input parameters;
% Example:
% [meanField,field,mat_file]=LESinverse(3,5,[],14,14,'gauss',1,0,1e-5,1e-2);

L=40;
if isempty(xy) && S==3 && R==5 && L==40
    load('direct 2D\xy.mat','xy_Optim');
    xy=xy_Optim;
elseif isempty(xy)
    xy=zeros(S+R,2); % first S rows contain x and y coordinates of sources;
    % the rest contain x and y coordinates of receivers
    %     linearCoord=4*2*L*rand(S+R,1);
    linearCoord=5+[4*L:(4*L)/S:8*L-4*L/S, 0:4*L/R:4*L-4*L/R]'; % location along one line with length=8*L
    %     linearCoord=[4*L:(4*L)/S:8*L-4*L/S, 4*L/R:(4*L-8*L/R)/(R-1):4*L-4*L/R]';
    % convert line coordinates to (x,y) coordinates
    ind=find(linearCoord<=L);
    xy(ind,2)=-L;
    xy(ind,1)=linearCoord(ind);
    ind=find(linearCoord>L & linearCoord<=3*L);
    xy(ind,1)=L;
    xy(ind,2)=linearCoord(ind)-L-L;
    ind=find(linearCoord>3*L & linearCoord<=5*L);
    xy(ind,2)=L;
    xy(ind,1)=-linearCoord(ind)+3*L+L;
    ind=find(linearCoord>5*L & linearCoord<=7*L);
    xy(ind,1)=-L;
    xy(ind,2)=-linearCoord(ind)+5*L+L;
    ind=find(linearCoord>7*L);
    xy(ind,2)=-L;
    xy(ind,1)=+linearCoord(ind)-7*L-L;
    xy1=xy;
    xy(1:S,:)=xy1(1:2:2*S,:);
    xy(S+1:S+R,:)=xy1([2:2:2*S 2*S+1:end],:);
end
Lx=[-L L];Ly=[-L L];
sigma.n=sigma_t;
sigma.x=sigma_x;
interp=1;
oldfile='No';
if exist('real_LES.mat','file')
    oldfile=questdlg('Previously saved file ''real_LES.mat'' is found. Would you like to use the saved data?',...
        '','Yes','No','Yes');
end
if isequal(oldfile,'No')
   % [tt,VX,VY,CC,T,c0,T0,vx0,vy0,xv,tau,xy]=real_LES(S,R,xy,sigman,mpl,atm,time,interp,showFig);
    [tt,VX,VY,CC,T,c0,T0,vx0,vy0,xv,tau,xy]=real_LES(S,R,xy,sigma,1,'u',[20:24],interp,0);
else
    load real_LES;
end
[Nray,Ntime]=size(tt);
yv=xv;
% ask for what frame to estimate (t0) and what frames to use for it (frames)
t0=1;
frames=1;
if Ntime>1
    disp(sprintf('The total number of frames is %d',Ntime)); 
    t0=input('What frame of fields to estimate?  ');
    if isempty(t0)
        disp(sprintf('The frame # %d will be estimated',Ntime)); 
        t0=Ntime;
    end
    frames=input('What frames of data to use?   ');
    if isempty(frames)
        disp('All frames of data will be used');
        frames=[1:Ntime];
    end
end
% prepare variables to launch showgraphs.m in order to see what are the
% original fields
mF.T=T0(t0);
mF.vx=vx0(t0);
mF.vy=vy0(t0);
mF.xy=xy;
mF.index=[1:S*R]';
mF.tt=tt(:,t0);
% mF.dtt=[];
T_original=T(:,:,t0);
vx_original=VX(:,:,t0);
vy_original=VY(:,:,t0);
for i=1:Ntime
    dVx(:,:,i)=VX(:,:,i)-vx0(i);
    dVy(:,:,i)=VY(:,:,i)-vy0(i);
    dT(:,:,i)=T(:,:,i)-T0(i);
    dCC(:,:,i)=CC(:,:,i)-c0(i);
end
sigma.vx=sqrt(mean(dVx(:).^2));
sigma.vy=sqrt(mean(dVy(:).^2));
sigma.T=sqrt(mean(dT(:).^2));
fluct.dT=dT(:,:,t0);
fluct.vx=dVx(:,:,t0);
fluct.vy=dVy(:,:,t0);
menuName=sprintf('Original fields at t=%g',t0);
showGraphs(menuName,S,R,mF,xv,xv,fluct,FigInterp,[],[]);
% ask whether the original fields are o.k. for the reconstruction; if no,
% the program will end.
buttonName=questdlg('Continue?');
if isequal(buttonName,'No')
    field=[];
    mat_file=[];
    meanField=[];
    return
end
% if the original fields are o.k. go ahead to reconstruct them
% mean fields reconstruction
i=0;
for j=[t0,frames]
    i=i+1;
    meanFields(i)=spatialMean(S,R,xy,tt(:,j),330,360,0);
%     [meanFields_T(i),meanFields_V(i)]=spatialMean_r(S,R,xy,tt(:,j),330,360,0);
end
U.x=[meanFields.vx]';
U.y=[meanFields.vy]';
U.sig=sqrt((sigma.vx^2+sigma.vy^2)/2);
% U_r.sig=U.sig;
% U_r.x=[meanFields_V.vx]';
% U_r.y=[meanFields_V.vy]';

c0_true=c0(t0);
T0_true=T0(t0)
vx0_true=vx0(t0)
vy0_true=vy0(t0)
meanField=meanFields(1)
% meanField_T=meanFields_T(1)
% meanField_V=meanFields_V(1)
if CnstrFlag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% constraints calculation
    grids=length(xv);
    % Laplacian for the T field
    k=0;
    gg1=zeros((grids-2)*(grids-2),grids*grids);
    for i=2:grids-1
        for j=2:grids-1
            pat=zeros(grids,grids);
            pat(i,j)=-4;
            pat([i-1,i+1],j)=1;
            pat(i,[j-1,j+1])=1;
            k=k+1;
            gg1(k,:)=pat(:)';
        end
    end
    % du/dx
    k=0;
    gg2=zeros((grids-2)*(grids-2),grids*grids);
    for i=2:grids-1
        for j=2:grids-1
            pat=zeros(grids,grids);
            pat(i-1,j)=-1;
            pat(i+1,j)=1;
            k=k+1;
            gg2(k,:)=pat(:)';
        end
    end
    % dv/dy
    k=0;
    gg3=zeros((grids-2)*(grids-2),grids*grids);
    for i=2:grids-1
        for j=2:grids-1
            pat=zeros(grids,grids);
            pat(i,j-1)=-1;
            pat(i,j+1)=1;
            k=k+1;
            gg3(k,:)=pat(:)';
        end
    end
    Cnstr.C=sparse([gg1 zeros((grids-2)*(grids-2),2*grids*grids);...
            zeros((grids-2)*(grids-2),grids*grids) gg2 gg3]);
    temp=[reshape(dT(:,:,t0),grids*grids,1);...
            reshape(dVx(:,:,t0),grids*grids,1);...
            reshape(dVy(:,:,t0),grids*grids,1)];
    Cnstr.B=sparse(Cnstr.C*temp);
else
    Cnstr=[];
end

% reconstruction of the fluctuations
EstFrameFlg='all';
SI=1;
[field,fields,mat_file]=t_stochastic_inverse2...
    ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    meanFields,U,tau,t0,frames,funType,EstFrameFlg,SI,Cnstr,'F',FigInterp,1);
save(mat_file,'mF','fluct','-append');

% calculation of actual and expected RMSE
T_rec=meanFields(1).T+field.dT;
vx_rec=meanFields(1).vx+field.vx;
vy_rec=meanFields(1).vy+field.vy;
res=T_rec(:)-T_original(:);
length_r=length(res);
T_actual_RMSE=sqrt(res.'*res/length_r)
T_expected_RMSE=sqrt(meanFields(1).std_dT^2+field.dTExpAverNMSE*sigma.T^2)
res=vx_rec(:)-vx_original(:);
vx_actual_RMSE=sqrt(res.'*res/length_r)
vx_expected_RMSE=sqrt(meanFields(1).std_dvx^2+field.vxExpAverNMSE*sigma.vx^2)
res=vy_rec(:)-vy_original(:);
vy_actual_RMSE=sqrt(res.'*res/length_r)
vy_expected_RMSE=sqrt(meanFields(1).std_dvy^2+field.vyExpAverNMSE*sigma.vy^2)

% sigma_T=sigma;
% sigma_T.vx=0;sigma_T.vy=0;
% % U_r=U;
% [field_T,fields_T]=t_stochastic_inverse2...
%     ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma_T,lT,lv,...
%     meanFields_T,U_r,tau,t0,frames,funType,EstFrameFlg,SI,Cnstr,'T',FigInterp,0);
% sigma_V=sigma;
% sigma_V.T=0;
% [field_V,fields_V]=t_stochastic_inverse2...
%     ('n',[],'n',[],S,R,Lx,Ly,xv,yv,sigma_V,lT,lv,...
%     meanFields_V,U_r,tau,t0,frames,funType,EstFrameFlg,SI,Cnstr,'V',FigInterp,0);
% 

%     menuName=sprintf('Unconstrained fields at t=%d',t0);
%     showGraphs(menuName,S,R,meanFields(1),xv,yv,field,FigInterp,[],[]);
%     menuName=sprintf('Unconstrained T fields at t=%d',t0);
%     showGraphs(menuName,S,R,meanFields_T(1),xv,yv,field_T,FigInterp,[],[]);
%     menuName=sprintf('Unconstrained V fields at t=%d',t0);
%     showGraphs(menuName,S,R,meanFields_V(1),xv,yv,field_V,FigInterp,[],[]);

