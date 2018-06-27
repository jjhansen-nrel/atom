function [field,fields,mat_file,R_dd0, R_md, G_1]=t_stochastic_inverse2...
    (ddType,R_dd0,mdType,R_md,S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
    meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig)
% time-dependent, constrained, and ordinary stochastic inversions for 2D locally frozen
% turbulence.
%
% Syntax:
% [field,fields,mat_file,R_dd0, R_md, G_1]=t_stochastic_inverse2...
%    (ddType,R_dd0,mdType,R_md,S,R,Lx,Ly,xv,yv,sigma,lT,lv,...
%    meanFields,U,tau,t0,frames,funType,estFrame,SI,Cnstr,InvType,interp,showFig);
%
% Inputs:
% ddType - an indicator of what Rdd0 matrix to use, a string; if ddType='o',
%   an old noise free Rdd0 matrix is used (this matrix should be specified in R_dd0;
%   if R_dd0=[], a menu to launch a ".mat" file where R_dd0 is stored will appear);
%   if ddType='n', a new matrix will be calculated;
% R_dd0 - an old noise free R_dd matrix; to make the program to use this matrix,
%   ddType must be equal 'o'; 
% mdType - an indicator of what Rmd matrix to use, a string; if mdType='o',
%   an old Rmd matrix is used (this matrix should be specified in R_md;
%   if R_md=[], a menu to launch a ".mat" file where R_md is stored will appear);
%   if mdType='n', a new matrix will be calculated;
% R_md - an old R_md matrix; to make the program to use this matrix,
%   mdType must be equal 'o';
% S - number of sound sources, integer;
% R - number of receivers, integer;
% Lx - a size of a tomographic area along x axis, e.g. Lx=[-40,40];
% Ly - a size of a tomographic area along y axis, e.g. Ly=[-40,40];
% xv - a vector of spatial points along x axis where the fields are recostructed;
% yv - a vector of spatial points along y axis where the fields are reconstructed;
% sigma - a struct variable with fields:
%   sigma.T - a STD of T fluctuations;
%   sigma.vx - a STD of vx fluctuations;
%   sigma.vy - a STD of vy fluctuations;
%   sigma.n - a STD of errors in travel time measurements, seconds;
%   sigma.x - a STD of errors in transducer positions, meters;
% lT, lv - correlation lengths of T and v fluctuations, meters;
% meanFields - a struct array of size [1,time+1]; the first element meanFields(1)
%   corresponds to the time at which the fields should be reconstructed;
%   the rest elements (2:time+1) correspond to moments of data; the fields
%   of meanFields are:
%   meanFields.T - the mean value of T field;
%   meanFields.c - the mean value of the Laplacian sound speed;
%   meanFields.vx - the mean value of vx field;
%   meanFields.vy - the mean value of vy field;
%   meanFields.xy - the x and y coordinated of sources and receivers; a
%   matrix [S+R, 2];
%   meanFields.data - a matrix of valid data and their indeces [data index]
%       for the reconstruction of the fluctuations;
%       data=-c_0^2*dtt, where dtt - travel times due to fluctations.
%   meanFields.std_dc - the standard deviation of the errors in c_0;
%   meanFields.std_dvx - the standard deviation of the errors in vx_0;
%   meanFields.std_dvy - the standard deviation of the errors in vy_0;
% U - a struct variable with fields:
%   U.x - the x component of mean wind velocity, m/s;
%   U.y - the y component of mean wind velocity, m/s;
%   U.sig - the standard deviation used in the formula for locally frozen
%   turbulence covariance function; for rigidly frozen turbulence U.sig=0.
% tau - time interval in seconds between two consecutive scans.  
% t0 - a time frame at which you would like to make a reconstruction, an integer;
% frames - a vector of frames of travel times which you would like to use as data;
%   for example, you have 5 times which correspond to 5 different frames of fields
%   and 5 vectors of travel times, then t0 can be 3 while frames can be [2,3,4].
% funType - a string which specifies what covariance functions should be
%   used in calculations; the options are:
%   'gauss' - Gaussian covariance functions;
%   'exp' - exponent;
%   'exp23' - exp[-(R/L)^(2/3)];
%   'vonk' - von Karman;
%   The calculation is the fastest with 'gauss' type;
% estFrame - a string which specifies what frames of fields should be
%   reconstructed. The options are:
%   'single' - fields only at time t0 will be reconstructed;
%   'all' - fields at times t0 and frames will be reconstructed;
% SI - a flag for ordinary stochastic inversion: SI=1 - ordinary SI will be
%   calculated; SI=0 - ordinary SI will not be calculated;
% Cnstr - constraints on the fluctuations, a struct variable with fields:
%   Cnstr.C - a matrix (usually, sparse matrix) of constraints: [number of
%   constraints, number of models]; the number of constraints should be
%   smaller than the number of models;
%   Cnstr.B - right hand-side values for the constraints, a vector: C*m=B;
% InvType - a flag what inversion to make; if InvType='T' - reconstruction
% of T fluctuations in the reciprocal transmission experiments
% (sigma.vx and sigma.vy must be 0); if InvType='V' - reconstruction of vx and vy
% fluctuations in the reciprocal transmission experiments (sigma.T must be
% 0); otherwise, both fields will be reconstructed.
% interp - interpolation of the displayed fields; output fields will not be interpolated;
%       interp = 0 - no interpolation; interp = 1 - one extra point will be added in the middle
%       of two adjacent points in space and time, ant etc.
% showFig - a flag to display graphics; if showFig = 1, graphics will be
%       displayed; if showFig = 0, the graphics will not be displayed;
%
% Outputs:
% field - the reconstructed fields at time t0, a struct variable with fields:
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
% fields - the reconstructed fields at times frames, a struct variable;
%   the fields of this vector are the same as of field, however, all matrices
%   are three dimensional: (x,y,t); fields = [] (empty) if estFrame = 'single';
% mat_file - a name of '.mat' file with the results; this file is
%   atomatically created in a sub-directory 'TDSI_Res' in the current dirctory
%   (if there is no such a sub-directory, it will be created); the name for the
%   first run of t_stochastic_inverse2 creates file '1.mat', the second -
%   '2.mat', and etc. After the 20th run, the file names will start wih 1
%   again, so that in the direcoty 'TDSI_Res' will be stored the results of
%   20 last runs of t_stochastic_inverse2; in this mat_file all vaiables
%   are stored, including constrained solution 'field_Cnstr' and 'fields_Cnstr',
%   and ordinary stochastic inversion 'field_SI' and 'fields_SI'; these variables
%   are struct variables with the same fields as 'field' and 'fields'; for each 
%   '.mat' file, a '.txt' file is also created with the same name containing a 
%   short description of the input parameters;
% R_dd0 - a noise free Rdd covariance matrix;
% R_md - R_md covariance matrix; if estFrame='all', R_md is a 3D array:
%   R_md(:,:,1) is R_md matrix for t0; R_md(:,:,2) is R_md matrix for the
%   first frame of data and so on.
% G_1 - a matrix G_1=R_md*pinv(R_dd), so that models m=G_1*data.

Ntime=length(frames);
if length(U.x)==1
    U.x=U.x*ones(Ntime+1,1);
end
if length(U.y)==1
    U.y=U.y*ones(Ntime+1,1);
end
sigmaT2=sigma.T^2;
sigmavx2=sigma.vx^2;
sigmavy2=sigma.vy^2;
sigmav2=sqrt(sigmavx2*sigmavy2);
sigman2=sigma.n^2;
sigmax2=sigma.x^2;
% sigmaT2=mean([meanFields(2:end).sig_dT].^2);
% sigmavx2=mean([meanFields(2:end).sig_dvx].^2);
% sigmavy2=mean([meanFields(2:end).sig_dvy].^2);
% sigmav2=sqrt(sigmavx2*sigmavy2);
% sigman2=mean([meanFields(2:end).std_tt].^2);
% sigma.T=sqrt(sigmaT2);
% sigma.vx=sqrt(sigmavx2);
% sigma.vy=sqrt(sigmavy2);
% sigma.n=sqrt(sigman2);
% U.x=[meanFields.vx]';U.y=[meanFields.vy]';
% U.sig=sqrt((sigmavx2+sigmavy2)/2);
% U.tau0=0;

if isempty(Cnstr)
    Cnstr_flag=0;
     buttonName='No';
else
    Cnstr_flag=1;
    text_dlg=sprintf(['Whould you like to calculate NMSE for the constrained reconstruction?\n'...
        'It takes a while...']);
    buttonName=questdlg(text_dlg,'','Yes for all frames','Yes for one frame','No','No');
    pause(1);
end
c0est=[meanFields.c];
T0est=[meanFields.T];
t0_frames=[t0,frames];
xy=meanFields(1).xy;
if isequal(estFrame,'all')
    ltf=Ntime+1;
else
    ltf=1;
    fields=[];
end
for i=1:S
    for j=1:R
        b=xy(S+j,1)-xy(i,1); % difference of x between the j-th receiver and i-th source 
        a=xy(S+j,2)-xy(i,2); % difference of y between the j-th receiver and i-th source
        length_r=sqrt(a*a+b*b);
        if length_r==0
            cos_phi=NaN;
            sin_phi=NaN;
        else
            cos_phi=b/length_r;
            sin_phi=a/length_r;
        end
        path(R*(i-1)+j).length=length_r;
        path(R*(i-1)+j).x0=xy(i,1);
        path(R*(i-1)+j).y0=xy(i,2);
        path(R*(i-1)+j).cos=cos_phi;
        path(R*(i-1)+j).sin=sin_phi;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% query for saved matrices 
if isequal(lower(ddType(1)),'o') && isempty(R_dd0)
    [FileName,PathName] = uigetfile('*.mat','Load R_dd0');
    if FileName~=0
        Old=load([PathName,FileName],'R_dd0','R_md');
        R_dd0=Old.R_dd0;
        R_md=Old.R_md;
    end
end   
if isequal(lower(mdType(1)),'o') && isempty(R_md)
    [FileName,PathName] = uigetfile('*.mat','Load R_{md}');
    if FileName~=0
        Old=load([PathName,FileName],'R_md');
        R_md=Old.R_md;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data vector formation
d_vec=[];
for i=1:Ntime
    d_vec=[d_vec;meanFields(i+1).data];
    ldata(i)=size(d_vec,1);
end
l_Rdd0=ldata(end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% R_dd0 matrix calculation
if isequal(lower(ddType(1)),'n') || isempty(R_dd0)
    R_dd0=zeros(l_Rdd0);
    B_dd=[];
    h_bar=waitbar(0,'R_{dd} matrix calculation','name',['0%: ' funType]);
    for i=1:l_Rdd0
        xMax1=path(d_vec(i,2)).length;
        ind=find(ldata>=i);
        t1=frames(ind(1));
        ind1=ind(1);
        c_0(1)=c0est(ind1+1);
        T_0(1)=T0est(ind1+1);
        path0=path(d_vec(i,2));
        for j=i:l_Rdd0
            xMax2=path(d_vec(j,2)).length;
            path1=path(d_vec(j,2));
            ind=find(ldata>=j);
            t2=frames(ind(1));
            ind2=ind(1);
            c_0(2)=c0est(ind2+1);
            T_0(2)=T0est(ind2+1);
            interval=(t2-t1)*tau;
            if ind1==ind2
                shift_x=0;
                shift_y=0;
                mean_U=0;
            elseif ind1<ind2
                shift_x=diff(frames(ind1:ind2))*(U.x(ind1+1:ind2)+U.x(ind1+2:ind2+1))*tau/2;
                shift_y=diff(frames(ind1:ind2))*(U.y(ind1+1:ind2)+U.y(ind1+2:ind2+1))*tau/2;
                mean_U=sqrt((shift_x/interval)^2+(shift_y/interval)^2);
            end
            %             if mean_U~=0
            %                 U.tau0=lv/mean_U;
            %                 % L2=(U.sig*U.tau0)^2*(sqrt(pi)*(t2-t1)*tau/U.tau0*erf((t2-t1)*tau/U.tau0)+exp(-((t2-t1)*tau/U.tau0)^2)-1);
            %                 L2=(U.sig*interval)^2*(1-5/6*(interval/U.tau0)^2);
            %             else
            L2=(U.sig*interval)^2;
            %             end
            sT2_eff=sigmaT2/(1+2*L2/lT^2)^(3/2);
            lT_eff=sqrt(lT*lT+2*L2);
            svx2_eff=sigmavx2/(1+2*L2/lv^2)^(3/2);
            lv_eff=sqrt(lv*lv+2*L2);
            svy2_eff=sigmavy2/(1+2*L2/lv^2)^(3/2);
            if ~isequal(lower(funType),'radial')
                path1.x0=path1.x0-shift_x;
                path1.y0=path1.y0-shift_y;
            else
                path1.x0=path1.x0-xv(t2);
                path1.y0=path1.y0-yv(t2);
                path0.x0=path(d_vec(i,2)).x0-xv(t1);
                path0.y0=path(d_vec(i,2)).y0-yv(t1);
            end
%             if isempty(B_dd)
%                 tf=[];
%             else
%                 tf=find(ismember(B_dd(:,1:5),[d_vec(i,2) d_vec(j,2) shift_x shift_y L2],'rows'));
%             end
%             if isempty(tf)
%                 tmp(1:5)=[d_vec(i,2) d_vec(j,2) shift_x shift_y L2];
                if isequal(lower(funType),'gauss')
                    R_dd0(i,j)=quad0('gaussCorr2D_an',0,xMax2,1e-8,[],path0,path1,...
                        c_0,T_0,lT_eff,lv_eff,sT2_eff,svx2_eff,svy2_eff);
                elseif isequal(lower(funType),'exp')
                    R_dd0(i,j)=dblquad('expCorr2D',0,xMax1,0,xMax2,1e-8,@quadl0,...
                        path0,path1,c_0,T_0,lT_eff,lv_eff,sT2_eff,svx2_eff,svy2_eff);
                elseif isequal(lower(funType),'exp23')
                    R_dd0(i,j)=dblquad('exp23Corr2D',0,xMax1,0,xMax2,1e-8,@quadl0,...
                        path0,path1,c_0,T_0,lT_eff,lv_eff,sT2_eff,svx2_eff,svy2_eff);
                elseif isequal(lower(funType),'vonk')
                    R_dd0(i,j)=dblquad('vonKCorr2D',0,xMax1,0,xMax2,1e-8,@quadl0,...
                        path0,path1,c_0,T_0,2*pi/lT_eff,2*pi/lv_eff,sT2_eff,svx2_eff,svy2_eff);
                elseif isequal(lower(funType),'radial')
                    xsi=path0.x0;
                    ysi=path0.y0;
                    xri=xsi+path0.length*path0.cos;
                    yri=ysi+path0.length*path0.sin;
                    xsj=path1.x0;
                    ysj=path1.y0;
                    xrj=xsj+path1.length*path1.cos;
                    yrj=ysj+path1.length*path1.sin;
                    R_dd0(i,j)=0.25*((xri^2-xsi^2)*(xrj^2-xsj^2)+(xrj^2-xsj^2)*(yri^2-ysi^2)+...
                        (xri^2-xsi^2)*(yrj^2-ysj^2)+(yri^2-ysi^2)*(yrj^2-ysj^2));
                end
%                tmp(6)=R_dd0(i,j);
%                 B_dd=[B_dd;tmp];
%             else
%                 R_dd0(i,j)=B_dd(tf,6);
%             end
%             clear tmp;
            R_dd0(j,i)=R_dd0(i,j);           
            waitbar(((i-1)*l_Rdd0+j)/(l_Rdd0*l_Rdd0),h_bar);
            set(h_bar,'name',sprintf('%d%%: %s',round(((i-1)*l_Rdd0+j)/(l_Rdd0*l_Rdd0)*100),funType));
        end
    end
    close(h_bar)
    clear B_dd
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% take noise into account and get R_dd
R_dd=R_dd0;
switch upper(InvType)
    case 'T'
        for i=1:l_Rdd0
            ind=find(ldata>=i);
            ind1=ind(1);
            for j=i:l_Rdd0
                ind=find(ldata>=j);
                ind2=ind(1);
                R_nn=0;     
                if d_vec(i,2)==d_vec(j,2)
                    R_nn=R_nn+2*sigmax2*c0est(ind1+1)*c0est(ind2+1);
                end
                if ind1==ind2
                    sig_dc02=meanFields(ind1+1).std_dc^2;
                    Li1=path(d_vec(i,2)).length;
                    Li2=path(d_vec(j,2)).length;
                    R_nn=R_nn+Li1*Li2*sig_dc02;
                end           
                if d_vec(i,2)==d_vec(j,2) && ind1==ind2
                    R_nn=R_nn+0.5*sigman2*c0est(ind1+1)^4;
                end
                R_dd(i,j)=R_dd(i,j)+R_nn;
                R_dd(j,i)=R_dd(i,j);
            end
        end
    case 'V'
        for i=1:l_Rdd0
            ind=find(ldata>=i);
            ind1=ind(1);
            for j=i:l_Rdd0
                ind=find(ldata>=j);
                ind2=ind(1);
                R_nn=0;     
                if d_vec(i,2)==d_vec(j,2)
                    R_nn=R_nn+2*sigmax2*(meanFields(ind1+1).vx*meanFields(ind2+1).vx+meanFields(ind1+1).vy*meanFields(ind2+1).vy);
                end
                if ind1==ind2
                    sig_dvx2=meanFields(ind1+1).std_dvx^2;
                    sig_dvy2=meanFields(ind1+1).std_dvy^2;
                    sx1=path(d_vec(j,2)).cos;
                    sy1=path(d_vec(j,2)).sin;
                    sx2=path(d_vec(i,2)).cos;
                    sy2=path(d_vec(i,2)).sin;
                    Li2=path(d_vec(i,2)).length;
                    Li1=path(d_vec(j,2)).length;
                    R_nn=R_nn+Li1*Li2*(sx1*sx2*sig_dvx2+sy1*sy2*sig_dvy2);
                end           
                if d_vec(i,2)==d_vec(j,2) && ind1==ind2
                    R_nn=R_nn+0.5*sigman2*c0est(ind1+1)^4;
                end
                R_dd(i,j)=R_dd(i,j)+R_nn;
                R_dd(j,i)=R_dd(i,j);
            end
        end
    otherwise
        for i=1:l_Rdd0
            ind=find(ldata>=i);
            ind1=ind(1);
            for j=i:l_Rdd0
                ind=find(ldata>=j);
                ind2=ind(1);
                R_nn=0;     
                if d_vec(i,2)==d_vec(j,2)
                    sx1=path(d_vec(j,2)).cos;
                    sy1=path(d_vec(j,2)).sin;
                    R_nn=R_nn+2*sigmax2*(c0est(ind1+1)*c0est(ind2+1)-(meanFields(ind1+1).vx*c0est(ind2+1)+...
                        meanFields(ind2+1).vx*c0est(ind1+1))*sx1-(meanFields(ind1+1).vy*c0est(ind2+1)+...
                        meanFields(ind2+1).vy*c0est(ind1+1))*sy1);
                end
                if ind1==ind2
                    sig_dvx2=meanFields(ind1+1).std_dvx^2;
                    sig_dvy2=meanFields(ind1+1).std_dvy^2;
                    sig_dc02=meanFields(ind1+1).std_dc^2;
                    sx1=path(d_vec(j,2)).cos;
                    sy1=path(d_vec(j,2)).sin;
                    sx2=path(d_vec(i,2)).cos;
                    sy2=path(d_vec(i,2)).sin;
                    Li2=path(d_vec(i,2)).length;
                    Li1=path(d_vec(j,2)).length;
                    R_nn=R_nn+Li1*Li2*(sig_dc02+sx1*sx2*sig_dvx2+sy1*sy2*sig_dvy2);
                end           
                if d_vec(i,2)==d_vec(j,2) && ind1==ind2
                    R_nn=R_nn+sigman2*c0est(ind1+1)^4;
                end
                R_dd(i,j)=R_dd(i,j)+R_nn;
                R_dd(j,i)=R_dd(i,j);
            end
        end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% R_md calculation
lyv=length(yv);
lxv=length(xv);
if isequal(lower(mdType(1)),'n') || isempty(R_md)
    R_md=zeros(3*lxv*lyv,l_Rdd0,ltf);
    h_bar=waitbar(0,'R_{md} matrices calculation','name',['0%: ' funType]);
    for j=1:lyv
        y2=yv(j);
        for k=1:l_Rdd0
            path1=path(d_vec(k,2));
            Li=path1.length;
            sin_phi=path1.sin;
            cos_phi=path1.cos;
            ind=find(ldata>=k);
            t2=frames(ind(1));
            ind2=ind(1)+1;
            for l=1:ltf
                t1=t0_frames(l);
                interval=(t2-t1)*tau;
                if t1==t2
                    shift_x=0;
                    shift_y=0;
                    mean_U=0;
                elseif t1<t2
                    ind_t=find(frames>t1 & frames<t2);
                    vec_t=[t1 frames(ind_t) t2];
                    U_x=[U.x(l); U.x(ind_t+1); U.x(ind2)];
                    U_y=[U.y(l); U.y(ind_t+1); U.y(ind2)];
                    shift_x=diff(vec_t)*(U_x(1:end-1)+U_x(2:end))*tau/2;
                    shift_y=diff(vec_t)*(U_y(1:end-1)+U_y(2:end))*tau/2;
                    %mean_U=sqrt((shift_x/interval)^2+(shift_y/interval)^2);
                elseif t1>t2
                    ind_t=find(frames<t1 & frames>t2);
                    vec_t=[t2 frames(ind_t) t1];
                    U_x=[U.x(ind2); U.x(ind_t+1); U.x(l)];
                    U_y=[U.y(ind2); U.y(ind_t+1); U.y(l)];
                    shift_x=-diff(vec_t)*(U_x(1:end-1)+U_x(2:end))*tau/2;
                    shift_y=-diff(vec_t)*(U_y(1:end-1)+U_y(2:end))*tau/2;
                    %mean_U=sqrt((shift_x/interval)^2+(shift_y/interval)^2);
                end
                %                 if mean_U~=0
                %                     U.tau0=lv/mean_U;
                % %                     L2=(U.sig*U.tau0)^2*(sqrt(pi)*(t2-t1)*tau/U.tau0*erf((t2-t1)*tau/U.tau0)+exp(-((t2-t1)*tau/U.tau0)^2)-1);
                %                     L2=(U.sig*interval)^2*(1-5/6*(interval/U.tau0)^2);
                %                 else
                L2=(U.sig*interval)^2;
                %                 end
                sT2_eff=sigmaT2/(1+2*L2/lT^2)^(3/2);
                lT_eff=sqrt(lT*lT+2*L2);
                svx2_eff=sigmavx2/(1+2*L2/lv^2)^(3/2);
                lv_eff=sqrt(lv*lv+2*L2);
                svy2_eff=sigmavy2/(1+2*L2/lv^2)^(3/2);
                sv2_eff=sqrt(svx2_eff*svy2_eff);
                if ~isequal(lower(funType),'radial')
                    path1.x0=path(d_vec(k,2)).x0-shift_x;
                    path1.y0=path(d_vec(k,2)).y0-shift_y;
                else
                    path1.x0=path(d_vec(k,2)).x0-xv(t2);
                    path1.y0=path(d_vec(k,2)).y0-yv(t2);
                end
                i0=1:lxv;
                if isequal(lower(funType),'gauss')
                    x2=xv(:);
                    i=1:lxv;
                    a=cos_phi*(path1.x0-x2)+sin_phi*(path1.y0-y2);
                    b=sin_phi*(path1.x0-x2)-cos_phi*(path1.y0-y2);
                    I1=0.5*sqrt(pi)*sT2_eff*lT_eff* exp(-b.*b/(lT_eff*lT_eff)).*(erf((Li+a)/lT_eff)-erf(a/lT_eff));
                    my=path1.y0-y2-a*sin_phi;
                    mx=path1.x0-x2-a*cos_phi;
                    e1=exp(-(Li+a)/lv_eff.*(Li+a)/lv_eff);
                    e2=exp(-a.*a/lv_eff/lv_eff);
                    e3=exp(-b.*b/(lv_eff*lv_eff));
                    dee3=e3.*(e1-e2);
                    bracee3=e3.*(-(Li+a)/(2*lv_eff).*e1+0.25*sqrt(pi)*erf((Li+a)/lv_eff)+...
                        a/(2*lv_eff).*e2-0.25*sqrt(pi)*erf(a/lv_eff));
                    I1v_core=0.5*sqrt(pi)*lv_eff*e3.* (erf((Li+a)/lv_eff)-erf(a/lv_eff));
                    I1vx=svx2_eff*I1v_core;
                    I1vy=I1v_core*svy2_eff;
                    I1v=I1v_core*sv2_eff;
                    c=a.*b*(cos_phi*cos_phi-sin_phi*sin_phi)-a.*a*sin_phi*cos_phi+(path1.x0-x2).*(path1.y0-y2);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%% I2
                    J1=I1vx;
                    J22=-my.*my/(lv_eff*lv_eff).*J1;
                    J21=-sin_phi*sin_phi*svx2_eff*lv_eff*bracee3;
                    J23=sin_phi*svx2_eff*dee3.*my;
                    J2=J21+J22+J23;
                    I2=J1+J2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I3
                    Y1=sin_phi*cos_phi*sv2_eff*lv_eff*bracee3;
                    Y2=0.5*sv2_eff*(cos_phi*cos_phi-sin_phi*sin_phi)*dee3.*b;
                    Y3=c/(lv_eff*lv_eff).*I1v;
                    I3=Y1+Y2+Y3;
                    I3_1=I3;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I4
                    Z1=I1vy;
                    Z22=-mx.*mx/(lv_eff*lv_eff).*Z1;
                    Z21=-cos_phi*cos_phi*svy2_eff*lv_eff*bracee3;
                    Z23=cos_phi*svy2_eff*dee3.*mx;
                    Z2=Z21+Z22+Z23;
                    I4=Z1+Z2;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                elseif isequal(lower(funType),'exp')
                    for i=1:lxv
                        I1(i,1)=sT2_eff*quadl0('exp_T_md_2D',0,Li,1e-8,[],path1,lT_eff,xv(i),yv(j));
                        I2(i,1)=svx2_eff*quadl0('exp_11_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                        I3(i,1)=sv2_eff*quadl0('exp_12_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                        I4(i,1)=svy2_eff*quadl0('exp_22_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                    end
                    I3_1=I3;
                elseif isequal(lower(funType),'exp23')
                    for i=1:lxv
                        I1(i,1)=sT2_eff*quadl0('exp23_T_md_2D',0,Li,1e-8,[],path1,lT_eff,xv(i),yv(j));
                        I2(i,1)=svx2_eff*quadl0('exp23_11_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                        I3(i,1)=sv2_eff*quadl0('exp23_12_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                        I4(i,1)=svy2_eff*quadl0('exp23_22_md_2D',0,Li,1e-8,[],path1,lv_eff,xv(i),yv(j));
                    end
                    I3_1=I3;
                elseif isequal(lower(funType),'vonk')
                    for i=1:lxv
                        I1(i,1)=sT2_eff*quadl0('vonK_T_md_2D',0,Li,1e-8,[],path1,2*pi/lT_eff,xv(i),yv(j));
                        I2(i,1)=svx2_eff*quadl0('vonK_11_md_2D',0,Li,1e-8,[],path1,2*pi/lv_eff,xv(i),yv(j));
                        I3(i,1)=sv2_eff*quadl0('vonK_12_md_2D',0,Li,1e-8,[],path1,2*pi/lv_eff,xv(i),yv(j));
                        I4(i,1)=svy2_eff*quadl0('vonK_22_md_2D',0,Li,1e-8,[],path1,2*pi/lv_eff,xv(i),yv(j));
                    end   
                    I3_1=I3;
                elseif isequal(lower(funType),'radial')
                    I1=zeros(lxv,1);
                    x2=xv(:)-xv(t1);
                    y2=yv(j)-yv(t1);
                    i=1:lxv;
                    I2(i,1)=0.5*x2*Li*(2*path1.x0+Li*cos_phi);
                    I3(i,1)=0.5*x2*Li*(2*path1.y0+Li*sin_phi);
                    I3_1(i,1)=repmat(0.5*y2*Li*(2*path1.x0+Li*cos_phi),lxv,1);
                    I4(i,1)=repmat(0.5*y2*Li*(2*path1.y0+Li*sin_phi),lxv,1);
                        
                end
                R_md(lxv*(j-1)+i0,k,l)=c0est(ind(1)+1)/(2*T0est(ind(1)+1))*I1;
                R_md(lxv*(j-1)+i0+lxv*lyv,k,l)=cos_phi*I2+sin_phi*I3;
                R_md(lxv*(j-1)+i0+2*lxv*lyv,k,l)=cos_phi*I3_1+sin_phi*I4;
                
                waitbar(((j-1)*ltf*l_Rdd0+(k-1)*ltf+l)/(lyv*ltf*l_Rdd0),h_bar);
                set(h_bar,'name',sprintf('%d%%: %s',round(...
                    ((j-1)*ltf*l_Rdd0+(k-1)*ltf+l)/(lyv*ltf*l_Rdd0)*100),funType));
            end % for l
        end % for k
        %             end
    end % for j
    close(h_bar);
end % if
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% models calculation
J=lxv*lyv;
SingValG=svd(R_dd);
temp=zeros(J,1);
SingValC=[];
for l=1:ltf
    G_1=R_md(:,:,l)*pinv(R_dd);
    m=G_1*d_vec(:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extraction of fields from the models    
    T=reshape(m(1:J),lxv,lyv);
    vx=reshape(m(J+1:2*J),lxv,lyv);
    vy=reshape(m(2*J+1:3*J),lxv,lyv);
    dC=c0est(l)/(2*T0est(l))*T;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expected NMSE calculation    
    covM=zeros(3*lxv,lyv);
    if sigmaT2~=0
        for i=1:J
            temp(i)=1-G_1(i,:)*R_md(i,:,l)'/sigmaT2;
        end
        covM(1:lxv,1:lyv)=reshape(temp,lxv,lyv);
    else
        covM(1:lxv,1:lyv)=zeros(lxv,lyv);
    end
    if sigmavx2~=0
        for i=1:J
            temp(i)=1-G_1(i+J,:)*R_md(i+J,:,l)'/sigmavx2;
        end
        covM(lxv+1:2*lxv,1:lyv)=reshape(temp,lxv,lyv);
    else
        covM(lxv+1:2*lxv,1:lyv)=zeros(lxv,lyv);
    end
    if sigmavy2~=0
        for i=1:J
            temp(i)=1-G_1(i+2*J,:)*R_md(i+2*J,:,l)'/sigmavy2;
        end
        covM(2*lxv+1:3*lxv,1:lyv)=reshape(temp,lxv,lyv);
    else
        covM(2*lxv+1:3*lxv,1:lyv)=zeros(lxv,lyv);
    end
    %%%%%%%%%%%%%%%%%%%%% creation of field (at t=t0) and fields (at t=frames) structures    
    if l==1
        field.dT=T;
        field.vx=vx;
        field.vy=vy;
        field.dC=dC;
        field.dTExpNMSE=covM(1:lxv,1:lyv);
        field.dTExpAverNMSE=mean(mean(field.dTExpNMSE));
        field.dTExpSTD=sqrt(field.dTExpNMSE*sigmaT2);
        field.dTExpAverSTD=mean(mean(field.dTExpSTD));
        field.vxExpNMSE=covM(lxv+1:2*lxv,1:lyv);
        field.vxExpAverNMSE=mean(mean(field.vxExpNMSE));
        field.vxExpSTD=sqrt(field.vxExpNMSE*sigmavx2);
        field.vxExpAverSTD=mean(mean(field.vxExpSTD));
        field.vyExpNMSE=covM(2*lxv+1:3*lxv,1:lyv);
        field.vyExpAverNMSE=mean(mean(field.vyExpNMSE));
        field.vyExpSTD=sqrt(field.vyExpNMSE*sigmavy2);
        field.vyExpAverSTD=mean(mean(field.vyExpSTD));
    else
        fields.dT(:,:,l-1)=T;
        fields.vx(:,:,l-1)=vx;
        fields.vy(:,:,l-1)=vy;
        fields.dC(:,:,l-1)=dC;
        fields.dTExpNMSE(:,:,l-1)=covM(1:lxv,1:lyv);
        fields.dTExpAverNMSE(l-1)=mean(mean(fields.dTExpNMSE(:,:,l-1)));
        fields.dTExpSTD(:,:,l-1)=sqrt(fields.dTExpNMSE(:,:,l-1)*sigmaT2);
        fields.dTExpAverSTD(l-1)=mean(mean(fields.dTExpSTD(:,:,l-1)));
        fields.vxExpNMSE(:,:,l-1)=covM(lxv+1:2*lxv,1:lyv);
        fields.vxExpAverNMSE(l-1)=mean(mean(fields.vxExpNMSE(:,:,l-1)));
        fields.vxExpSTD(:,:,l-1)=sqrt(fields.vxExpNMSE(:,:,l-1)*sigmavx2);
        fields.vxExpAverSTD(l-1)=mean(mean(fields.vxExpSTD(:,:,l-1)));
        fields.vyExpNMSE(:,:,l-1)=covM(2*lxv+1:3*lxv,1:lyv);
        fields.vyExpAverNMSE(l-1)=mean(mean(fields.vyExpNMSE(:,:,l-1)));
        fields.vyExpSTD(:,:,l-1)=sqrt(fields.vyExpNMSE(:,:,l-1)*sigmavy2);
        fields.vyExpAverSTD(l-1)=mean(mean(fields.vyExpSTD(:,:,l-1)));
    end    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% saving results
if exist([cd '\TDSI_Res'],'dir')
    files=dir([cd '\TDSI_Res\*.mat']);
    if ~isempty(files)
        [mi,i]=max(datenum({files(:).date}));
        curr_num=str2num(files(i).name(1:end-4)); %#ok<ST2NM>
        name=[cd '\TDSI_Res\' num2str(mod(curr_num,20)+1) '.txt'];
    else
        name=[cd '\TDSI_Res\1.txt'];
    end
else
    mkdir('TDSI_Res');
    name=[cd '\TDSI_Res\1.txt'];
end
mat_file=[name(1:end-4) '.mat'];
save(mat_file,'R_dd','R_dd0','R_md','field','fields','S','R','Lx','Ly','xv','yv','sigma','lT','lv',...
    'xy','U','tau','t0','frames','funType','estFrame','Cnstr','Cnstr_flag','SI','d_vec','ldata','c0est','T0est',...
    'SingValG','SingValC','lxv','lyv','meanFields','interp','ltf','J');
fid=fopen(name,'wt');
fprintf(fid,[' S=%d \n R=%d \n Lx=[%s] m \n Ly=[%s] m \n x_step=%g m \n y_step=%g m \n x_length=%g pt \n '...
        'y_length=%g pt \n sigma_noise=%g s \n '...
        'sigma_T=%g K \n sigma_vx=%g m/s \n sigma_vy=%g m/s \n sigma_xy=%g m \n lT=%g m \n lv=%g m \n Ux=[%s] m/s \n '...
        'Uy=[%s] m/s \n sigma_U=%g m/s \n tau=%g s \n t0=%d \n frames=[%s] \n rays=[%s] \n '...
        'funType=''%s'' \n estFrame=''%s'' \n Constraints=%d \n SI=%d'],S,R,sprintf('%g ',Lx),sprintf('%g ',Ly),...
    xv(2)-xv(1),yv(2)-yv(1),lxv,lyv,sigma.n,sigma.T,sigma.vx,sigma.vy,sigma.x,lT,lv,sprintf('%g ',U.x'),sprintf('%g ',U.y'),...
    U.sig,tau,t0,sprintf('%g ',frames),sprintf('%g ',diff([0 ldata])),funType,estFrame,Cnstr_flag,SI);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   graphics
if showFig
    menuName=sprintf('Unconstrained fields at t=%d',t0);
    %(menuName,S,R,meanFields,xv,yv,fluct,SingValG,SingValC)
    showGraphs(menuName,S,R,meanFields(1),xv,yv,field,interp,SingValG,SingValC);
    if ~isempty(fields)
        menuName='Unconstrained fields at specified frames';
        showGraphs(menuName,S,R,meanFields(2:end),xv,yv,fields,interp,SingValG,SingValC);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stochastic inversion
if SI
    pause(1);
    if ~(length(frames)==1 && t0==frames)
        save ('mat_file_name.mat','mat_file','showFig','Cnstr_flag','buttonName');
        clear all;
        load ('mat_file_name.mat');
        load(mat_file,'ldata','d_vec','R_dd','R_md','c0est','T0est','lxv','lyv','sigma','J','ltf','t0','frames');
        sigmaT2=sigma.T*sigma.T;
        sigmavx2=sigma.vx*sigma.vx;
        sigmavy2=sigma.vy*sigma.vy;
        l_d=[0 ldata];
        fields_SI=[];
        for l=1:ltf
            flg=1;
            i0=l-1;
            if l==1
                i0=find(frames==t0);
                if isempty(i0)
                    disp(sprintf('Cannot calculate SI for %g frame',t0));
                    flg=0;
                    field_SI=[];
                end
            end
            if flg
                Rdd_SI=R_dd(l_d(i0)+1:l_d(i0+1),l_d(i0)+1:l_d(i0+1));
                Rmd_SI=R_md(:,l_d(i0)+1:l_d(i0+1),l);
                d_SI=d_vec(l_d(i0)+1:l_d(i0+1),1);
                G1_SI=Rmd_SI*pinv(Rdd_SI);
                m=G1_SI*d_SI;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% extraction of fields from the models    
                T=reshape(m(1:J),lxv,lyv);
                vx=reshape(m(J+1:2*J),lxv,lyv);
                vy=reshape(m(2*J+1:3*J),lxv,lyv);
                dC=c0est(l)/(2*T0est(l))*T;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% expected NMSE calculation    
                covM=zeros(3*lxv,lyv);
                if sigmaT2~=0
                    for i=1:J
                        temp(i)=1-G1_SI(i,:)*Rmd_SI(i,:)'/sigmaT2;
                    end
                    covM(1:lxv,1:lyv)=reshape(temp,lxv,lyv);
                else
                    covM(1:lxv,1:lyv)=zeros(lxv,lyv);
                end
                if sigmavx2~=0
                    for i=1:J
                        temp(i)=1-G1_SI(i+J,:)*Rmd_SI(i+J,:)'/sigmavx2;
                    end
                    covM(lxv+1:2*lxv,1:lyv)=reshape(temp,lxv,lyv);
                else
                    covM(lxv+1:2*lxv,1:lyv)=zeros(lxv,lyv);
                end
                if sigmavy2~=0
                    for i=1:J
                        temp(i)=1-G1_SI(i+2*J,:)*Rmd_SI(i+2*J,:)'/sigmavy2;
                    end
                    covM(2*lxv+1:3*lxv,1:lyv)=reshape(temp,lxv,lyv);
                else
                    covM(2*lxv+1:3*lxv,1:lyv)=zeros(lxv,lyv);
                end
                %%%%%%%%%%%%%%%%%%%%% creation of field (at t=t0) and fields (at t=frames) structures    
                if l==1
                    field_SI.dT=T;
                    field_SI.vx=vx;
                    field_SI.vy=vy;
                    field_SI.dC=dC;
                    field_SI.dTExpNMSE=covM(1:lxv,1:lyv);
                    field_SI.dTExpAverNMSE=mean(mean(field_SI.dTExpNMSE));
                    field_SI.dTExpSTD=sqrt(field_SI.dTExpNMSE*sigmaT2);
                    field_SI.dTExpAverSTD=mean(mean(field_SI.dTExpSTD));
                    field_SI.vxExpNMSE=covM(lxv+1:2*lxv,1:lyv);
                    field_SI.vxExpAverNMSE=mean(mean(field_SI.vxExpNMSE));
                    field_SI.vxExpSTD=sqrt(field_SI.vxExpNMSE*sigmavx2);
                    field_SI.vxExpAverSTD=mean(mean(field_SI.vxExpSTD));
                    field_SI.vyExpNMSE=covM(2*lxv+1:3*lxv,1:lyv);
                    field_SI.vyExpAverNMSE=mean(mean(field_SI.vyExpNMSE));
                    field_SI.vyExpSTD=sqrt(field_SI.vyExpNMSE*sigmavy2);
                    field_SI.vyExpAverSTD=mean(mean(field_SI.vyExpSTD));
                else
                    fields_SI.dT(:,:,l-1)=T;
                    fields_SI.vx(:,:,l-1)=vx;
                    fields_SI.vy(:,:,l-1)=vy;
                    fields_SI.dC(:,:,l-1)=dC;
                    fields_SI.dTExpNMSE(:,:,l-1)=covM(1:lxv,1:lyv);
                    fields_SI.dTExpAverNMSE(l-1)=mean(mean(fields_SI.dTExpNMSE(:,:,l-1)));
                    fields_SI.dTExpSTD(:,:,l-1)=sqrt(fields_SI.dTExpNMSE(:,:,l-1)*sigmaT2);
                    fields_SI.dTExpAverSTD(l-1)=mean(mean(fields_SI.dTExpSTD(:,:,l-1)));
                    fields_SI.vxExpNMSE(:,:,l-1)=covM(lxv+1:2*lxv,1:lyv);
                    fields_SI.vxExpAverNMSE(l-1)=mean(mean(fields_SI.vxExpNMSE(:,:,l-1)));
                    fields_SI.vxExpSTD(:,:,l-1)=sqrt(fields_SI.vxExpNMSE(:,:,l-1)*sigmavx2);
                    fields_SI.vxExpAverSTD(l-1)=mean(mean(fields_SI.vxExpSTD(:,:,l-1)));
                    fields_SI.vyExpNMSE(:,:,l-1)=covM(2*lxv+1:3*lxv,1:lyv);
                    fields_SI.vyExpAverNMSE(l-1)=mean(mean(fields_SI.vyExpNMSE(:,:,l-1)));
                    fields_SI.vyExpSTD(:,:,l-1)=sqrt(fields_SI.vyExpNMSE(:,:,l-1)*sigmavy2);
                    fields_SI.vyExpAverSTD(l-1)=mean(mean(fields_SI.vyExpSTD(:,:,l-1)));
                end % if l==1    
            end % if flg
        end %for l=1:ltf
        save(mat_file,'field_SI','fields_SI','-append');
        clear all  
        load('mat_file_name.mat');
        if showFig
            load(mat_file,'S','R','Lx','Ly','xv','yv','meanFields',...
                'SingValG','SingValC','frames','interp','t0','field_SI','fields_SI');
            menuName=sprintf('Unconstrained SI fields at t=%d',t0);
            %(menuName,S,R,meanFields,xv,yv,fluct,SingValG,SingValC)
            showGraphs(menuName,S,R,meanFields(1),xv,yv,field_SI,interp,SingValG,SingValC);
            if ~isempty(fields_SI)
                menuName='Unconstrained SI fields at specified frames';
                showGraphs(menuName,S,R,meanFields(2:end),xv,yv,fields_SI,interp,SingValG,SingValC);
            end
        end
    end % if ~(length(frames)==1 & t0==frames)
end % if SI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% constrained solution
if Cnstr_flag
    pause(1);
    save ('mat_file_name.mat','mat_file','showFig','buttonName');
    clear all;
    load ('mat_file_name.mat');
    load(mat_file,'Cnstr');
    C=Cnstr.C1;
    if isempty(C) || all(all(C==0))
        load(mat_file,'field','fields');
        field_Cnstr.dT=field.dT;
        field_Cnstr.dC=field.dC;
        field_Cnstr.dTExpNMSE=field.dTExpNMSE;
        field_Cnstr.dTExpAverNMSE=field.dTExpAverNMSE;
        field_Cnstr.dTExpSTD=field.dTExpSTD;
        field_Cnstr.dTExpAverSTD=field.dTExpAverSTD;
        for i=1:length(fields)
            fields_Cnstr(i).dT=fields(i).dT;
            fields_Cnstr(i).dC=fields(i).dC;
            fields_Cnstr(i).dTExpNMSE=fields(i).dTExpNMSE;
            fields_Cnstr(i).dTExpAverNMSE=fields(i).dTExpAverNMSE;
            fields_Cnstr(i).dTExpSTD=fields(i).dTExpSTD;
            fields_Cnstr(i).dTExpAverSTD=fields(i).dTExpAverSTD;
        end
        
    else
        load(mat_file,'sigma','c0est','T0est','J','lxv','lyv','ltf');
        B=Cnstr.B1;
        SingValC.T =[];
        CT=C';
        sigmaT2=sigma.T*sigma.T;
        fields_Cnstr=[];
        load(mat_file,'R_md','R_dd','d_vec');
        p_R_dd=pinv(R_dd);
        clear R_dd;
        G_1=p_R_dd*d_vec(:,1);
        for l=1:ltf
%             G_1=R_md(:,:,l)*p_R_dd;
            m=R_md(1:J,:,l)*G_1;
            %      m=m-C_core*(C*m-B);
            if isempty(B) || all(B==0)
                m_T=m-CT*((C*CT)\(C*m));
            else
                m_T=m-CT*((C*CT)\(C*m-B));
            end
            if isequal(buttonName,'Yes for all frames')
                temp_T=NMSE_Cnstr_2D_1(p_R_dd,R_md(:,:,l),J,C,B,sigmaT2,0,0); 
            elseif isequal(buttonName,'No')
                temp_T=repmat(0,J,1);
            else
                if l==1
                    temp_T=NMSE_Cnstr_2D_1(p_R_dd,R_md(:,:,l),J,C,B,sigmaT2,0,0);
                    
                else
                    temp_T=repmat(0,J,1);
                    
                end
            end
            clear G_1;
            if l==1
                field_Cnstr.dT=reshape(m_T,lxv,lyv);
                field_Cnstr.dC=c0est(l)/(2*T0est(l))*field_Cnstr.dT;
                field_Cnstr.dTExpNMSE=reshape(temp_T,lxv,lyv);
                field_Cnstr.dTExpAverNMSE=mean(mean(field_Cnstr.dTExpNMSE));
                field_Cnstr.dTExpSTD=sqrt(field_Cnstr.dTExpNMSE*sigmaT2);
                field_Cnstr.dTExpAverSTD=mean(mean(field_Cnstr.dTExpSTD));
            else
                fields_Cnstr(l-1).dT=reshape(m_T,lxv,lyv);
                fields_Cnstr(l-1).dC=c0est(l)/(2*T0est(l))*fields_Cnstr(l-1).dT;
                fields_Cnstr(l-1).dTExpNMSE=reshape(temp_T,lxv,lyv);
                fields_Cnstr(l-1).dTExpAverNMSE=mean(mean(fields_Cnstr(l-1).dTExpNMSE));
                fields_Cnstr(l-1).dTExpSTD=sqrt(fields_Cnstr(l-1).dTExpNMSE*sigmaT2);
                fields_Cnstr(l-1).dTExpAverSTD=mean(mean(fields_Cnstr(l-1).dTExpSTD));
            end    % if l==1
        end % for l=1:ltf
    end % if isempty(C) | all(C==0)
    clear C;
    C=Cnstr.C2;
    if isempty(C) || all(all(C==0))
        load(mat_file,'field','fields');
        field_Cnstr.vx=field.vx;
        field_Cnstr.vy=field.vy;
        field_Cnstr.vxExpNMSE=field.vxExpNMSE;
        field_Cnstr.vyExpNMSE=field.vyExpNMSE;
        field_Cnstr.vxExpAverNMSE=field.vxExpAverNMSE;
        field_Cnstr.vxExpSTD=field.vxExpSTD;
        field_Cnstr.vxExpAverSTD=field.vxExpAverSTD;
        field_Cnstr.vyExpAverNMSE=field.vyExpAverNMSE;
        field_Cnstr.vyExpSTD=field.vyExpSTD;
        field_Cnstr.vyExpAverSTD=field.vyExpAverSTD;
        for i=1:length(fields)
            fields_Cnstr(i).vx=fields(i).vx;
            fields_Cnstr(i).vy=fields(i).vy;
            fields_Cnstr(i).vxExpNMSE=fields(i).vxExpNMSE;
            fields_Cnstr(i).vyExpNMSE=fields(i).vyExpNMSE;
            fields_Cnstr(i).vxExpAverNMSE=fields(i).vxExpAverNMSE;
            fields_Cnstr(i).vxExpSTD=fields(i).vxExpSTD;
            fields_Cnstr(i).vxExpAverSTD=fields(i).vxExpAverSTD;
            fields_Cnstr(i).vyExpAverNMSE=fields(i).vyExpAverNMSE;
            fields_Cnstr(i).vyExpSTD=fields(i).vyExpSTD;
            fields_Cnstr(i).vyExpAverSTD=fields(i).vyExpAverSTD;
        end
        
    else
        load(mat_file,'sigma','c0est','T0est','J','lxv','lyv','ltf');
        B=Cnstr.B2;
        SingValC.V=[]; 
        CT=C';
        sigmavx2=sigma.vx*sigma.vx;
        sigmavy2=sigma.vy*sigma.vy;
        fields_Cnstr=[];
        if ~exist('R_md','var') || ~exist('p_R_dd','var') || ~exist('d_vec','var')
            load(mat_file,'R_md','R_dd','d_vec');
            p_R_dd=pinv(R_dd);
            clear R_dd;
        end
        G_1=p_R_dd*d_vec(:,1);
        for l=1:ltf
%             G_1=R_md(:,:,l)*p_R_dd;
            m=R_md(:,:,l)*G_1;
            %      m=m-C_core*(C*m-B);
            if isempty(B) || all(B==0)
                m_V=m(J+1:3*J)-CT*((C*CT)\(C*m(J+1:3*J)));
            else
                m_V=m(J+1:3*J)-CT*((C*CT)\(C*m(J+1:3*J)-B));
            end
            if isequal(buttonName,'Yes for all frames')
                [temp_T,temp_vx,temp_vy]=NMSE_Cnstr_2D_1(p_R_dd,R_md(:,:,l),J,C,B,0,sigmavx2,sigmavy2); 
                clear temp_T;
            elseif isequal(buttonName,'No')
                temp_vx=repmat(0,J,1);
                temp_vy=temp_vx;
            else
                if l==1
                    [temp_T,temp_vx,temp_vy]=NMSE_Cnstr_2D_1(p_R_dd,R_md(:,:,l),J,C,B,0,sigmavx2,sigmavy2);
                    clear temp_T;
                else
                    temp_vx=repmat(0,J,1);
                    temp_vy=temp_vx;
                end
            end
            clear G_1;
            if l==1
                field_Cnstr.vx=reshape(m_V(1:J),lxv,lyv);
                field_Cnstr.vy=reshape(m_V(J+1:2*J),lxv,lyv);
                field_Cnstr.vxExpNMSE=reshape(temp_vx,lxv,lyv);
                field_Cnstr.vyExpNMSE=reshape(temp_vy,lxv,lyv);
                field_Cnstr.vxExpAverNMSE=mean(mean(field_Cnstr.vxExpNMSE));
                field_Cnstr.vxExpSTD=sqrt(field_Cnstr.vxExpNMSE*sigmavx2);
                field_Cnstr.vxExpAverSTD=mean(mean(field_Cnstr.vxExpSTD));
                field_Cnstr.vyExpAverNMSE=mean(mean(field_Cnstr.vyExpNMSE));
                field_Cnstr.vyExpSTD=sqrt(field_Cnstr.vyExpNMSE*sigmavy2);
                field_Cnstr.vyExpAverSTD=mean(mean(field_Cnstr.vyExpSTD));
            else
                fields_Cnstr(l-1).vx=reshape(m_V(1:J),lxv,lyv);
                fields_Cnstr(l-1).vy=reshape(m_V(J+1:2*J),lxv,lyv);
                fields_Cnstr(l-1).vxExpNMSE=reshape(temp_vx,lxv,lyv);
                fields_Cnstr(l-1).vyExpNMSE=reshape(temp_vy,lxv,lyv);
                fields_Cnstr(l-1).vxExpAverNMSE=mean(mean(fields_Cnstr(l-1).vxExpNMSE));
                fields_Cnstr(l-1).vxExpSTD=sqrt(fields_Cnstr(l-1).vxExpNMSE*sigmavx2);
                fields_Cnstr(l-1).vxExpAverSTD=mean(mean(fields_Cnstr(l-1).vxExpSTD));
                fields_Cnstr(l-1).vyExpAverNMSE=mean(mean(fields_Cnstr(l-1).vyExpNMSE));
                fields_Cnstr(l-1).vyExpSTD=sqrt(fields_Cnstr(l-1).vyExpNMSE*sigmavy2);
                fields_Cnstr(l-1).vyExpAverSTD=mean(mean(fields_Cnstr(l-1).vyExpSTD));
            end    % if l==1
        end % for l=1:ltf
    end % if isempty(C) | all(C==0)
    
    save(mat_file,'field_Cnstr','fields_Cnstr','SingValC','-append');
    clear all  
    load('mat_file_name.mat');
    if showFig
        load(mat_file,'S','R','Lx','Ly','xv','yv','meanFields','field_Cnstr','fields_Cnstr',...
            'SingValG','SingValC','frames','interp','t0');
        if length(frames)==1 && frames==t0
            menuName=sprintf('Constrained SI at t=%d',t0);
        else
            menuName=sprintf('Constrained TDSI at t=%d',t0);
        end
        %(menuName,S,R,meanFields,xv,yv,fluct,SingValG,SingValC)
        showGraphs(menuName,S,R,meanFields(1),xv,yv,field_Cnstr,interp,SingValG,SingValC);
        if ~isempty(fields_Cnstr)
            menuName='Constrained TDSI at specified frames';
            showGraphs(menuName,S,R,meanFields(2:end),xv,yv,fields_Cnstr,interp,SingValG,SingValC);
        end
    end
end % if Cnstr_flag
if nargout>0 && ~exist('field','var')
    switch nargout
        case 1
            load(mat_file,'field')
        case {2,3}
            load(mat_file,'field','fields');
        case 4
            load(mat_file,'field','fields','R_dd0');
        case 5
            load(mat_file,'field','fields','R_dd0', 'R_md');
        case 6
            load(mat_file,'field','fields','R_dd0', 'R_md', 'G_1');
    end
end