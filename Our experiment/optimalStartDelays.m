function [delays,fVal,s]=optimalStartDelays
maxDel(1)=-20; % in ms
maxDel(2)=20;
step=4;
Nsig=8;
fVal=zeros(1,Nsig);
% delays0=zeros(Nsig,1); % in ms
 delays0=[120;100;200;0;160;200;40;140]; % fVal=15.432; i=4; fval(second min)=17.9 ms
 % delays0=[128;108;208;0;164;200;40;140]; % fVal = 11 ms; i=4; fval(second min)=18.5 ms
% delays0=[124; 104; 204; 0; 160; 200; 40; 144]; % fVal= 17.311; i=4;
%   delays0=[180;0;160;0;180;160;40;40]; % fVal=15.4; i=2;
% delays0=[180;0;156;0;180;160;36;36]; % fVal=15.424; i=2;
% delays0=[3.7082
%         16.44
%        10.063
%        5.0636
%        13.796
%        15.919
%        8.2183];
%    delay0=[11.98
%        23.638
%         20.32
%        10.938
%        12.589
%        25.264
%        20.563];
% delays0=[10.156
%     5.4725
%     13.863
%     8.3049
%     6.3704
%     10.641
%     16.44];
% delays0=[14.077
%     9.3931
%     17.784
%     10.494
%     18.166
%     19.745
%     19.177];
% delays0=[21.389
%     5.8093
%     15.568
%     13.648
%     27.151
%     26.893
%     24.433];
% delay0=[31.071
%     9.6257
%     14.673
%     22.799
%     21.95
%     37.47
%     32.68];
% delays=zeros(Nsig-1,Nsig);% in ms
% LB=zeros(Nsig-1,1);
% UB=maxDel*ones(Nsig-1,1);
xy=[40.6579 79.7960
    17.5507 1.1259
    77.8724 12.7153
    57.1284 0
    0.8467  67.3008
    0       24.7025
    79.7780 55.0969
    73.0493 79.1313];
xy=xy-40;% the centered coordinates

% array of speakers
xyS=xy;
xyS(4,1)=xyS(4,1)+0.4;
xyS(5,2)=xyS(5,2)-0.4;
xyS(6,2)=xyS(6,2)-0.4;
xyS(7,2)=xyS(7,2)+0.4;
xyS(8,1)=xyS(8,1)-0.4;
% array of receivers
xyR=xy;
xyR(1,1)=xyR(1,1)+0.8;
xyR(2,1)=xyR(2,1)-0.8;
xyR(3,2)=xyR(3,2)-0.8;
xyR(4,1)=xyR(4,1)-0.4;
xyR(5,2)=xyR(5,2)+0.4;
xyR(6,2)=xyR(6,2)+0.4;
xyR(7,2)=xyR(7,2)-0.4;
xyR(8,1)=xyR(8,1)+0.4;
%% new tower's indexing
tmp=xyS;
newIndex=1:8;
oldIndex=indexMapN2O(newIndex);
xyS(newIndex,:)=tmp(oldIndex,:);
tmp=xyR;
xyR(newIndex,:)=tmp(oldIndex,:);
tt=zeros(Nsig,Nsig);
c0=343;
ell=zeros(Nsig,Nsig);
for i=1:Nsig
    for j=1:Nsig
        % i-th source j-th receiver
        ell(i,j)=sqrt((xyS(i,:)-xyR(j,:))*(xyS(i,:)-xyR(j,:)).');
        tt(i,j)=ell(i,j)/c0;
        %         delay(i,j)=round(tt(i,j)/dt);
    end
end
tt=1000*tt;% in ms

%% fmin
% options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
% for i=1:Nsig
%     [delays(:,i),fVal(i)] = fmincon(@(s) PS(s,i,Nsig,tt),delays0,[],[],[],[],LB,UB,[],options);
% end
% fVal=-fVal;
% [~,i]=max(fVal);
% s=delays(:,i);
% [~,md]=PS(s,i,Nsig,tt);
% md

%% brute force
delays=zeros(Nsig,1);
md0=zeros(Nsig,1);
for i=4:4%1:Nsig
    [delays(:,i),fVal(i),md0(:,i)]=bruteForce(maxDel,step,delays0,i,Nsig,tt);
end
[~,i]=max(fVal);
s=delays(:,i)
md=md0(:,i)


function [crit,md]=PS(s,i,Nsig,tt)
speakerDelayTime(1:i-1,1)=s(1:i-1);
speakerDelayTime(i,1)=0;
speakerDelayTime(i+1:Nsig,1)=s(i:Nsig-1);
ttAsynch=tt+repmat(speakerDelayTime,1,Nsig);
md=zeros(Nsig,1);
for i=1:Nsig
    st=sort(ttAsynch(:,i));% sorted tt for i-th mic over all speakers
    md(i)=min(diff(st));
end
% crit=-prod(md);
crit=-min(md);

function [speakerDelayTime0,crit0,md0]=bruteForce(delMax,step,delays0,i,Nsig,tt)

speakerDelayTime=zeros(Nsig,1);
speakerDelayTime0=zeros(Nsig,1);
crit0=0;
speakerDelayTime(i)=0;
% delMax=delMax+1;
md=zeros(Nsig,1);
% N=fix(delMax/step)+1;
range=delMax(1):step:delMax(2);
for i1=range
    for i2=range
        for i3=range
            for i4=range
                for i5=range
                    for i6=range
                        for i7=range
%                             s=[i1;i2;i3;i4;i5;i6;i7]-1;
%                             s=s*step;
                            s=[delays0(1:i-1);delays0(i+1:end)];
                            s=s+[i1;i2;i3;i4;i5;i6;i7];
                            speakerDelayTime(1:i-1)=s(1:i-1);
                            speakerDelayTime(i+1:Nsig)=s(i:end);
                            ttAsynch=tt+repmat(speakerDelayTime,1,Nsig);
                            for j=1:Nsig
                                st=sort(ttAsynch(:,j));% sorted tt for i-th mic over all speakers
                                md(j)=min(diff(st));
                            end
%                              crit=prod(md);
                            [~,j]=min(md);
                            crit=min(md([1:j-1,j+1:end]));
                            if crit>crit0
                                crit0=crit
                                speakerDelayTime0=speakerDelayTime
                                md0=md
                            end
                        end
                    end
                end
            end
        end
    end
end

function [speakerDelayTime0,crit0,md0]=bruteForce1(delMax,i,Nsig,tt)

speakerDelayTime=zeros(Nsig,1);
speakerDelayTime0=zeros(Nsig,1);
crit0=0;
speakerDelayTime=[7; 9; 6; 0; 15; 22; 23; 9];
delMax=delMax+1;
md=zeros(Nsig,1);
for i1=1:1
    for i2=1:delMax
        for i3=1:delMax
            for i4=1:1
                for i5=1:delMax
                    for i6=1:delMax
                        for i7=1:1
                            s=[i1;i2;i3;i4;i5;i6;i7]-1;
                            speakerDelayTime(1:i-1)=s(1:i-1);
                            speakerDelayTime(i+1:Nsig)=s(i:end);
                            ttAsynch=tt+repmat(speakerDelayTime,1,Nsig);
                            for j=1:Nsig
                                st=sort(ttAsynch(:,j));% sorted tt for i-th mic over all speakers
                                md(j)=min(diff(st));
                            end
                            % crit=prod(md);
                            crit=min(md);
                            if crit>crit0
                                crit0=crit
                                speakerDelayTime0=speakerDelayTime
                                md0=md
                            end
                        end
                    end
                end
            end
        end
    end
end

function oldIndex=indexMapN2O(newIndex)
newOrder=[1,8,7,3,4,2,6,5,9,10];
oldIndex=newOrder(newIndex);