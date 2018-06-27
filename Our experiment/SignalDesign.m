function [y,ttAsynch]=SignalDesign(sigType)
% SignalDesign i sa simulation tool for the BAO acoustic tomography array.
% sigType={'chirp' or 'tone'} constructs either chirp (wideband) or tone
% (narrow band) signals, both windowed by the hanning window.
% y is a matrix of [Ntime,Nsignals] containing signals (Nsignals=8) o flength ntime samples.

% speakerDelayTime=[7; 9; 6; 0; 15; 22; 23; 9]*1e-3;% in new indexing, in ms
% speakerDelayTime=[20.87; 0; 1; 0; 5.5; 25.24; 0; 0.32]*1e-3;
% speakerDelayTime=[0; 13; 8; 0; 0; 9; 12; 3]*1e-3;

%% This is an optimal delays; criterion: product of minimal time differences on all mics
% speakerDelayTime=[12; 24; 20; 11; 13; 0; 25; 21]*1e-3;
%% This is an optimal delays; criterion: min of minimal time differences on mics.
% speakerDelayTime=[31; 10; 15; 23; 22; 37; 0; 33]*1e-3;
%% brute force maxD=10 ms; crit: min
speakerDelayTime=[0; 10; 6; 0; 0; 6; 10; 0]*1e-3;
%% 100 ms
 speakerDelayTime=[80; 60; 60; 30; 0; 30; 90; 80]*1e-3;% 11.265

 %%
fs=100e3; % sampling frequency (Hz)
dt=1/fs; % discretization time
Nsig=8;
t_max=5.8e-3; % duration of the impulse (s)
% dt=0.5*0.25*1e-4; % discretization time (s);
% fs=1/dt; % sampling frequency (Hz)
t=(0:dt:t_max).'; % time grid
lt=length(t);
y=zeros(lt,Nsig);
% assign colors for speakers
col=[0 0 1
    1 0 0
    0 1 0
    1 0 1
    1 0.8 0
    0 0 0
    0 0.85 0.85
    0.85 0.33 0.1];
close all
pause(1)
switch sigType
    case 'chirp'
        F1=500;
        %         F2=7500;
        F2=2500;
        
        y(:,1)=chirp(t,F1,t_max,F2,'li',90);
        % show the chirp signal
        figure;
        plot(t*1000,y(:,1))
        axis tight
        xlabel('t (ms)','fontsize',12,'fontweight','bold');
        ylabel('Chirp signal','fontsize',12,'fontweight','bold');
        title(sprintf('T=%3.1f (ms) Fs=%3.1f (kHz) F1=%3.1f (kHz) F2=%3.1f (kHz)',...
            t_max*1000,fs*1e-3,F1*1e-3,F2*1e-3),'fontsize',12,'fontweight','bold');
        % make pivot functions
        %% sinusoids
        y1=zeros(lt,Nsig);
        fPivot=linspace(F1,F2,Nsig);
        for n=1:Nsig
            %     y1(:,n)=sin((n+3)*2*pi/t_max*t);
            y1(:,n)=sin(2*pi*fPivot(n)*t);
        end
        %% step-like functions
        % y1=zeros(lt,Nsig-1);
        % for n=1:Nsig-1
        %     l=lt/2^n;
        %     for m=1:2^n
        %         y1((m-1)*l+1:m*l,n)=(-1)^m;
        %     end
        % end
        %% show the pivot functions
%         figure;
%         for n=1:Nsig
%             subplot(2,Nsig/2,n)
%             plot(t*1000,y1(:,n))
%             axis tight
%             xlabel('t (ms)','fontsize',12,'fontweight','bold');
%             ylabel(['Signal ',num2str(n)],'fontsize',12,'fontweight','bold');
%             if n==1
%                 title('Pivot functions','fontsize',12,'fontweight','bold');
%             end
%         end
        
        %% construct the final signals
%         y=0.5*(y(:,1)*ones(1,Nsig)+y1);
        y=y(:,1)*ones(1,Nsig);
    case 'tone'
        F1=1000;%Hz
        F2=5000;
        sF1=indexMapO2N([1 2 3 5]);
        sF2=indexMapO2N([4 6 7 8]);
        y(:,sF1)=sin(2*pi*F1*t)*ones(1,length(sF1));
        y(:,sF2)=sin(2*pi*F2*t)*ones(1,length(sF2));
end

w=hanning(lt);
y=y.*(w*ones(1,Nsig));
% show final signals
figure;
for n=1:Nsig
    subplot(2,Nsig/2,n)
    plot(t*1000,y(:,n),'color',col(n,:))
    axis tight
    xlabel('t (ms)','fontsize',12,'fontweight','bold');
    ylabel(['Signal ',num2str(n)],'fontsize',12,'fontweight','bold');
    if n==1
        title('Windowed signals: hanning','fontsize',12,'fontweight','bold');
    end
end

%% correlations
% C=xcorr(y,lt/2,'coeff');
% figure;imagesc(1:Nsig^2,((1:length(C))-lt/2-1)*dt*1e3,C);colorbar
% xlabel('Signal pair','fontsize',12,'fontweight','bold');
% ylabel('Delay time (ms)','fontsize',12,'fontweight','bold');
% for n=1:Nsig
%     figure;
%     plot(((1:length(C))-lt/2-1)*dt*1e3,C(:,(n-1)*Nsig+1:n*Nsig))
%     xlabel('Delay time (ms)','fontsize',12,'fontweight','bold');
%     ylabel('Normilized cross-correlation','fontsize',12,'fontweight','bold');
%     for j=1:Nsig
%         txt{j}=sprintf('%d - %d',n,j);
%     end
%     legend(txt,0);
% end
% typical travel times
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
xyS(newIndex,:)=tmp(oldIndex,:) %#ok<NOPRT>
% xyS(2,:)=tmp(8,:);
% xyS(3,:)=tmp(7,:);
% xyS(4,:)=tmp(3,:);
% xyS(5,:)=tmp(4,:);
% xyS(6,:)=tmp(2,:);
% xyS(7,:)=tmp(6,:);
% xyS(8,:)=tmp(5,:);
% recivers
tmp=xyR;
xyR(newIndex,:)=tmp(oldIndex,:) %#ok<NOPRT>
% xyR(2,:)=tmp(8,:);
% xyR(3,:)=tmp(7,:);
% xyR(4,:)=tmp(3,:);
% xyR(5,:)=tmp(4,:);
% xyR(6,:)=tmp(2,:);
% xyR(7,:)=tmp(6,:);
% xyR(8,:)=tmp(5,:);
%%
tt=zeros(Nsig,Nsig);
delay=zeros(Nsig,Nsig);
c0=343;
ell=zeros(Nsig,Nsig);
for i=1:Nsig
    for j=1:Nsig
        % i-th source j-th receiver
        ell(i,j)=sqrt((xyS(i,:)-xyR(j,:))*(xyS(i,:)-xyR(j,:)).');
        tt(i,j)=ell(i,j)/c0;
        delay(i,j)=round(tt(i,j)/dt);
    end
end
% delay1=delay(indexMapO2N((1:Nsig)),:);
delay=delay+repmat(round(speakerDelayTime/dt),1,Nsig);
ttAsynch=tt+repmat(speakerDelayTime,1,Nsig);
% delay1=delay1(indexMapN2O((1:Nsig)),:);
% corrections for different starting times
% delay(indexMapO2N([1 5 7 8]),:)=delay(indexMapO2N([1 5 7 8]),:)+round(9e-3/dt);
% delay(indexMapO2N(2),:)=delay(indexMapO2N(2),:)+round(19e-3/dt);
% delay(indexMapO2N(4),:)=delay(indexMapO2N(4),:)+round(15e-3/dt);
% delay(indexMapO2N(6),:)=delay(indexMapO2N(6),:)+round(18e-3/dt);
% %% new corrections
% delay(indexMapO2N(1),:)=delay(indexMapO2N(1),:)-round(2e-3/dt);
% delay(indexMapO2N(2),:)=delay(indexMapO2N(2),:)+round(3e-3/dt);
% delay(indexMapO2N(6),:)=delay(indexMapO2N(6),:)+round(5e-3/dt);
% delay(indexMapO2N(7),:)=delay(indexMapO2N(7),:)-round(3e-3/dt);
% % corrections for asynchronic transmission
% ttAsynch=tt;
% ttAsynch(indexMapO2N([1 5 7 8]),:)=ttAsynch(indexMapO2N([1 5 7 8]),:)+9e-3;
% ttAsynch(indexMapO2N(2),:)=ttAsynch(indexMapO2N(2),:)+19e-3;
% ttAsynch(indexMapO2N(4),:)=ttAsynch(indexMapO2N(4),:)+15e-3;
% ttAsynch(indexMapO2N(6),:)=ttAsynch(indexMapO2N(6),:)+18e-3;
% %% new corrections
% ttAsynch(indexMapO2N(1),:)=ttAsynch(indexMapO2N(1),:)-2e-3;
% ttAsynch(indexMapO2N(2),:)=ttAsynch(indexMapO2N(2),:)+3e-3;
% ttAsynch(indexMapO2N(6),:)=ttAsynch(indexMapO2N(6),:)+5e-3;
% ttAsynch(indexMapO2N(7),:)=ttAsynch(indexMapO2N(7),:)-3e-3;
%% new tower's indexing
%%
% figure
% stem(ell(:))
% xlabel('Ray path','fontsize',12,'fontweight','bold');
% ylabel('Path length (m)','fontsize',12,'fontweight','bold');
% figure
% imagesc(ell)
% axis xy
% colorbar
% xlabel('Sources','fontsize',12,'fontweight','bold');
% ylabel('Receivers','fontsize',12,'fontweight','bold');
%% Signals on microphones
nrec=lt+max(max(delay));
tRec=(0:nrec-1)*dt*1000;
switch sigType
    case 'chirp'
        recS=zeros(nrec,Nsig);
        for i=1:Nsig
            figure
            plot(tRec,zeros(nrec,1),'b');
            hold on
            for j=1:Nsig
                %  recSTag{i,j}=sprintf('%d-%d',i,j);
                %                 recS(delay(j,i)+1:delay(j,i)+lt,i)=recS(delay(j,i)+1:delay(j,i)+lt,i)+ y(:,j);
                plot((delay(j,i)+1:delay(j,i)+lt)*dt*1000,y(:,j),'color',col(j,:));
            end
            set(gca,'Xlim',[tRec(1) tRec(end)],'XTick',0:10:tRec(end),'XGrid','on');
            xlabel('t (ms)','fontsize',12,'fontweight','bold');
            ylabel(['Microphone ',num2str(i)],'fontsize',12,'fontweight','bold');
            title('Recorded signals','fontsize',12,'fontweight','bold');
        end
        
        %         for n=1:Nsig
        %             figure
        %             plot(tRec,recS(:,n))
        %             xlabel('t (ms)','fontsize',12,'fontweight','bold');
        %             ylabel(['Microphone ',num2str(n)],'fontsize',12,'fontweight','bold');
        %             set(gca,'Xlim',[tRec(1) tRec(end)],'XTick',0:10:tRec(end),'XGrid','on');
        %             title('Recorded signals','fontsize',12,'fontweight','bold');
        %
        %         end
        
    case 'tone'
        recSF1=zeros(nrec,Nsig);
        for i=1:Nsig
            for j=sF1
                %  recSTag{i,j}=sprintf('%d-%d',i,j);
                recSF1(delay(j,i)+1:delay(j,i)+lt,i)=recSF1(delay(j,i)+1:delay(j,i)+lt,i)+ y(:,j);
            end
            %             figure
            %             plot(tRec,recSF1(:,i),'r')
            %             xlabel('t (ms)','fontsize',12,'fontweight','bold');
            %             ylabel(['Microphone ',num2str(i)],'fontsize',12,'fontweight','bold');
            %             axis tight
            %             title('Recorded signals','fontsize',12,'fontweight','bold');
            %             set(gca,'xtick',0:10:tRec(end),'xgrid','on')
        end
        recSF2=zeros(nrec,Nsig);
        for i=1:Nsig
            for j=sF2
                %  recSTag{i,j}=sprintf('%d-%d',i,j);
                recSF2(delay(j,i)+1:delay(j,i)+lt,i)=recSF2(delay(j,i)+1:delay(j,i)+lt,i)+ y(:,j);
            end
            %             figure
            %             plot(tRec,recSF2(:,i),'b')
            %             xlabel('t (ms)','fontsize',12,'fontweight','bold');
            %             ylabel(['Microphone ',num2str(i)],'fontsize',12,'fontweight','bold');
            %             axis tight
            %             title('Recorded signals','fontsize',12,'fontweight','bold');
            %             set(gca,'xtick',0:10:tRec(end),'xgrid','on')
        end
        for i=1:Nsig
            figure
            plot(tRec,recSF1(:,i),'r',tRec,recSF2(:,i),'b')
            xlabel('t (ms)','fontsize',12,'fontweight','bold');
            ylabel(['Microphone ',num2str(i)],'fontsize',12,'fontweight','bold');
            axis tight
            title('Recorded signals','fontsize',12,'fontweight','bold');
            set(gca,'xtick',0:10:tRec(end),'xgrid','on')
        end
end
% reconstruction of the travel times
% y1=[y;zeros(nrec-lt,Nsig)];
% ttR=zeros(Nsig,Nsig);

%%

for i=1:Nsig
    %     figure
    %     for j=1:Nsig
    %         [C,lags]=xcorr(recS(:,i),y1(:,j),'coeff');
    %         [mC,tmp]=max(C);
    %         ttR(j,i)=lags(tmp)*dt;
    %         subplot(2,Nsig/2,j)
    %         plot(lags*dt*1000,C)
    %         xlim([-2 2]+lags(tmp)*dt*1000)
    %         txt=sprintf('%d-%d. Max=%.2f at t=%.1f',i,j,mC,ttR(j,i)*1000);
    %         xlabel('\tau (ms)','fontsize',12,'fontweight','bold');
    %         ylabel('Correlation function','fontsize',12,'fontweight','bold');
    %         title(txt,'fontsize',12,'fontweight','bold');
    %     end
    figure
    %     plot(1:Nsig,tt(:,i)*1000,'*b',1:Nsig,ttR(:,i)*1000,'or')
    switch sigType
        case 'tone'
            for j=1:Nsig
                if j==indexMapO2N(1) || j==indexMapO2N(2) || j==indexMapO2N(3) || j==indexMapO2N(5)
                    stem(j,ttAsynch(j,i)*1000,'r')
                    hold on
                else
                    stem(j,ttAsynch(j,i)*1000,'b')
                    hold on
                end
            end
        case 'chirp'
            for j=1:Nsig
                    stem(j,ttAsynch(j,i)*1000,'color',col(j,:))
                    hold on
            end
    end
    set(gca,'ylim',[0 tRec(end)],'xlim',[0.5 Nsig+0.5],'ytick',0:10:tRec(end),'ygrid','on')
    xlabel('Speaker index','fontsize',12,'fontweight','bold');
    ylabel('Travel time and delay (ms)','fontsize',12,'fontweight','bold');
    title(['Microphone ',num2str(i)],'fontsize',12,'fontweight','bold')
    %     legend('Actual','Reconstructed');
end

function oldIndex=indexMapN2O(newIndex)
newOrder=[1,8,7,3,4,2,6,5,9,10];
oldIndex=newOrder(newIndex);

function newIndex=indexMapO2N(oldIndex)
oldOrder=[1,6,4,5,8,7,3,2,9,10];
newIndex=oldOrder(oldIndex);