function xy=arrayCalibration(t,xy0)
% Estimates (x,y) coordinates of effective point-sources to minimize the
% errors of reconstruction in the mean fields

tMean=mean(t,2);% temporary mean of travel times
S=3; % number of speakers
% there are S*2 variables - the effective source coordinates
% The initial guess is the measured coordinates
X0=xy0(1:S,:);
X0=X0(:);
xyR=xy0(S+1:end,:);
% define the function of mismatch
f=@(x) ttMismatch(x,tMean,xyR);
X=fminsearch(f,X0);
xy=[reshape(X,S,2);xyR];

function f=ttMismatch(x,tMean,xyR)
lx=length(x);
S=lx/2;
R=size(xyR,1);
xy=[reshape(x,S,2);xyR];
Lo=300;
Hi=400;
meanFields=spatialMean(S,R,xy,tMean,Lo,Hi,0);
f=meanFields.std_dc^2+meanFields.std_dvx^2+meanFields.std_dvy^2;
