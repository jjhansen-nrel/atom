function nmse = fluctErrors
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%  load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Our experiment\TDSI_Res\18.mat')
% load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Our experiment\TDSI_Res\8.mat')
%  load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Our experiment\2008 August Andreas\TDSI_Res\3.mat');
 load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Inverse 2D\TDSI_Res\3.mat')
%  load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Our experiment\2012_7_2 squared_off_Hanning\1438-1445\Travel times filterFlg 1 upSampleFlg 1\TDSI_Res\1.mat');
%  load('C:\Users\u4rrdsnv\Documents\MATLAB\work\Tomography\Our experiment\2008 July\2131-2133\Travel times filterFlg 1 upSampleFlg 1\TDSI_Res\8.mat')
path=createPath(S,R,xy);
l_Rdd0=size(R_dd0,1);
% mF=meanFields;% save original meanFields errors
%% no errors in the data
field=expectedErrors(lxv,lyv,sigma,R_md,R_dd0);
nmse.NoiseFreeData=[field.dTExpAverNMSE field.vxExpAverNMSE field.vyExpAverNMSE]
%% errors in travel times
sigma.x=0;
for i=1:length(meanFields)
    meanFields(i).std_dvx=0;
    meanFields(i).std_dvy=0;
    meanFields(i).std_dc=0;
end
%
N=50;
std_t=linspace(0,0.1,N);
std_t=std_t*1e-3; % in sec now
for i=1:N
    sigma.n=std_t(i);
    R_nn=calcRdd(meanFields,l_Rdd0,ldata,d_vec,sigma,path,c0est);
    R_dd=R_dd0+R_nn;
    field=expectedErrors(lxv,lyv,sigma,R_md,R_dd);
    nmse.tt(i,:)=[field.dTExpAverNMSE field.vxExpAverNMSE field.vyExpAverNMSE];
end
%
figure;
plot(std_t*1000,nmse.tt,'linewidth',2)
xlabel('Std in travel times (ms)','fontsize',12,'fontweight','bold');
ylabel('Average NMSE in fluctuations','fontsize',12,'fontweight','bold');
legend('\DeltaT','\Deltavx','\Deltavy',0)
grid
drawnow
%% errors in transducer coordinates
sigma.n=0;
%
std_x=linspace(0,5,N);
std_x=std_x*1e-2; % in (m) now
for i=1:N
    sigma.x=std_x(i);
    R_nn=calcRdd(meanFields,l_Rdd0,ldata,d_vec,sigma,path,c0est);
    R_dd=R_dd0+R_nn;
    field=expectedErrors(lxv,lyv,sigma,R_md,R_dd);
    nmse.dx(i,:)=[field.dTExpAverNMSE field.vxExpAverNMSE field.vyExpAverNMSE];
end
%
figure;
plot(std_x*100,nmse.dx,'linewidth',2)
xlabel('Std in coordinates (cm)','fontsize',12,'fontweight','bold');
ylabel('Average NMSE in fluctuations','fontsize',12,'fontweight','bold');
legend('\DeltaT','\Deltavx','\Deltavy',0)
grid
drawnow
%% errors in mean field estimations
sigma.n=0;
sigma.x=0;
%
std_v=linspace(0,0.5,N);

for j=1:N
    for i=1:length(meanFields)
        meanFields(i).std_dvx=std_v(j);
        meanFields(i).std_dvy=std_v(j);
        meanFields(i).std_dc=std_v(j);
    end
    R_nn=calcRdd(meanFields,l_Rdd0,ldata,d_vec,sigma,path,c0est);
    R_dd=R_dd0+R_nn;
    field=expectedErrors(lxv,lyv,sigma,R_md,R_dd);
    nmse.meanFields(j,:)=[field.dTExpAverNMSE field.vxExpAverNMSE field.vyExpAverNMSE];
end
figure;
plot(std_v,nmse.meanFields,'linewidth',2)
xlabel('Std in mean fields (m/s)','fontsize',12,'fontweight','bold');
ylabel('Average NMSE in fluctuations','fontsize',12,'fontweight','bold');
legend('\DeltaT','\Deltavx','\Deltavy',0)
grid
drawnow
%% errors in all factors
%
% std_v=linspace(0,0.5,N);
L0=mean([path(:).length]);

std_v=c0est(1)/L0*sqrt(2*std_x.^2+c0est(1)^2*std_t.^2);
for j=1:N
    for i=1:length(meanFields)
        meanFields(i).std_dvx=std_v(j);
        meanFields(i).std_dvy=std_v(j);
        meanFields(i).std_dc=std_v(j);
    end
    sigma.n=std_t(j);
    sigma.x=std_x(j);
    R_nn=calcRdd(meanFields,l_Rdd0,ldata,d_vec,sigma,path,c0est);
    R_dd=R_dd0+R_nn;
    field=expectedErrors(lxv,lyv,sigma,R_md,R_dd);
    nmse.all(j,:)=[field.dTExpAverNMSE field.vxExpAverNMSE field.vyExpAverNMSE];
end
figure;
plot(nmse.all,'linewidth',2)
xlabel('Index','fontsize',12,'fontweight','bold');
ylabel('Average NMSE in fluctuations','fontsize',12,'fontweight','bold');
legend('\DeltaT','\Deltavx','\Deltavy',0)
grid
drawnow

function field=expectedErrors(lxv,lyv,sigma,R_md,R_dd)
l=1;
J=lxv*lyv;
G_1=R_md(:,:,l)*pinv(R_dd);
sigmaT2=sigma.T^2;
sigmavx2=sigma.vx^2;
sigmavy2=sigma.vy^2;
% sigmav2=sqrt(sigmavx2*sigmavy2);
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


function R_dd=calcRdd(meanFields,l_Rdd0,ldata,d_vec,sigma,path,c0est)
sigman2=sigma.n^2;
sigmax2=sigma.x^2;
R_dd=zeros(l_Rdd0);
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

function path=createPath(S,R,xy)
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
