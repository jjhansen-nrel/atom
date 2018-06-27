function [meanField,delays]=spatialMean_serial_t(S,R,xy,tt,Lo,Hi,showFig)
% spatialMean_serial_t works like spatialMean_serial, that is, uses several
% time scans to estimate the delays and the mean fields, but the latters
% dpend on time.
% The input tt is a matrix of
% size [S*R, time] rather than a vector [S*R,1]. All travel times measured
% at different times are used for the reconstruction of mean T(t) and V(t)
% fields. It assumes that all data are collected at the same physical time
% The output meanField is a struct vector of length
% time rather than a struct variable.

[~, Ntime]=size(tt);
d=[];
SR=S*R;
G=zeros(SR*Ntime,3*Ntime+SR);
index=zeros(Ntime+1,1);
Li=cell(Ntime,1);t=Li;in=Li;
for i=1:Ntime
    [Li{i},t{i},in{i},theta]=ttFiltr(S,R,xy,tt(:,i),Lo,Hi,showFig);
    d=[d;t{i}./Li{i}];
    l=length(Li{i});
%     G(l0+1:l0+l,(i-1)*3+1)=1;
%     G(l0+1:l0+l,(i-1)*3+[1:2])=-theta;
    G((i-1)*SR+in{i},(i-1)*3+1)=1;
    G((i-1)*SR+in{i},(i-1)*3+[2:3])=-theta;
        for j=1:l
            G((i-1)*SR+in{i}(j),in{i}(j)+3*Ntime)=1/Li{i}(j);
        end
    index(i+1)=length(d);
end
lLi=index(Ntime+1);
%% clear zero rows
tmp=[];
for i=1:SR*Ntime
    if all(G(i,:)==0)
        tmp=[tmp;i];
    end
end
if ~isempty(tmp)
    G(tmp,:)=[];
end
%%
G_1=pinv(G);
aver=G_1*d;
delays=aver(3*Ntime+1:end);

er=d-G*aver;
% er1=d-G*av;
sig2=er'*er/(lLi-size(G,2));
std_a=sqrt(sig2*diag(G_1*G_1'));


for i=1:Ntime
    dtt=Li{i}.*er(index(i)+1:index(i+1));
    c0est=1/aver((i-1)*3+1);
    meanField(i).c=c0est;
    %meanField.dc=dc;
    c02=c0est*c0est;
    T0est=c02/343/343*293;
    meanField(i).T=T0est;
    %meanField.dT=dT;
    vx0est=aver((i-1)*3+2)*c02;
    meanField(i).vx=vx0est;
    %meanField.dvx=dvx;
    vy0est=aver((i-1)*3+3)*c02;
    meanField(i).vy=vy0est;
    %meanField.dvy=dvy;
    dc=c02*std_a((i-1)*3+1);
    meanField(i).std_dc=dc;
    dc2=2*c0est*dc;
    dT=dc2/343/343*293;
    meanField(i).std_dT=dT;
    dvx=c02*std_a((i-1)*3+2);%+aver(2)*dc2;
    meanField(i).std_dvx=dvx;
    dvy=c02*std_a((i-1)*3+3);%+aver(3)*dc2;
    meanField(i).std_dvy=dvy;
    %%%%%%%%%%%%%%%%%%%%%%%%% misc
    % meanField.std_tt=sig_tt;
    % meanField.std_x=sig_x;
    meanField(i).dtt=dtt;
    meanField(i).tt=t{i};
    meanField(i).index=in{i};
    meanField(i).data=[-c0est*c0est*dtt in{i}];
    meanField(i).xy=xy;
    meanField(i).Li=Li{i};
end