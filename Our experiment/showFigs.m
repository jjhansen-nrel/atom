function [tt_out,tt]=showFigs(tt,s_r,S,M,dts,filterFlg,upSampleFlg,filename,outlierFlg)
% displays travel times and corresponding signals on receivers. One click 
% above or below a marker indicating travel time will open a separate
% figure with the corresponding signal on the receiver. Do not click ON the
% marker; it will not display the signal.
% The program can be called without any arguments. Then, a file with the extracted
% travel times should be specified (see tt_extractor.m). If called with
% arguments, then the inputs and outputs are as follows.
%
% Inputs:
% tt - travel times (ms), a matrix of size [Nrays, M*N_minutes], where M indicates
% the number of scans in 1 min. For consecutive scanning during 1 min, M=120 (a 0.5 s
% interval between scans).
% s_r - a structure containing the signals on receivers.
% S - the number of sources.
% M - number of scans in 1 min.
% dts - the time interval between two samples in the digital signal; it is
% inverse to the sampling frequency.
% filterFlg  - a filtering flag; 1-the filtration was "on" during the travel
% time extraxtion; 0 - the filtration was "off".
% upSampleFlg - a upsampling flag; 1-upsampling was "on" during the travel
% time extraction; 0 - "off".
% filename - the name of the file being processed.
% outlierFlg - a flag turning "on"/"off" the automatic outlier remover;
% curently, there is no satisfactory automatic procedure, so, set this 0.
% 
% Outputs:
% tt_out - a matrix of the travel times in the format [Nrays, M/S];
% tt - original mastrix of travel times.

if nargin==0
    [file_name,path_name]=uigetfile('*.mat','Pick a file with tt');
    if file_name==0
        return
    end
    load([path_name,file_name]);
    dts=dt; %#ok<NODEF>
    filename=[file_name(4:7),'/',file_name(8:9),'/',file_name(10:11),' ',file_name(12:13),':',file_name(14:15),':',file_name(16:17)];
    outlierFlg=0;
end
dt=60/M; % record interval in sec.
Choice=1;
filter=num2str(filterFlg);
resam=num2str(upSampleFlg);
U.N=S;
tt_out=zeros(size(tt,1),M/S);
thresh=1.4;
% K=3;
while Choice~=16
    Choice=menu(filename,'S1R1','S1R2','S1R3','S1R4', 'S1R5',...
        'S2R1','S2R2','S2R3','S2R4', 'S2R5',...
        'S3R1','S3R2','S3R3','S3R4', 'S3R5','Exit');
    switch Choice
        case 1 % S1R1
            figure;
            plot((0:S:M-1)*dt,tt(1,1:S:M),'o');
            if outlierFlg
                %tt_out(1,:)=outlier_remover(tt(1,1:S:M),0.15);
                  tt_out(1,:)=outlier_remover3(tt(1,1:S:M),thresh);
%                  tt_out(1,:)=outlier_remover2(tt(1,1:S:M),K,thresh);
                hold on
                plot((0:S:M-1)*dt,tt_out(1,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,1:S:M,1));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S1R1: ','Filtration=',filter,' Resampling=',resam];
            U.S=1;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 2 % S1R2
            figure;
            plot((0:S:M-1)*dt,tt(2,1:S:M),'o');
            if outlierFlg
%                 tt_out(2,:)=outlier_remover(tt(2,1:S:M),0.15);
                tt_out(2,:)=outlier_remover3(tt(2,1:S:M),thresh);
%                 tt_out(2,:)=outlier_remover2(tt(2,1:S:M),K,thresh);
                hold on
                plot((0:S:M-1)*dt,tt_out(2,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,1:S:M,2));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S1R2: ','Filtration=',filter,' Resampling=',resam];
            U.S=1;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 3 % S1R3
            figure;
            plot((0:S:M-1)*dt,tt(3,1:S:M),'o');
            if outlierFlg
%                 tt_out(3,:)=outlier_remover(tt(3,1:S:M),0.15);
                tt_out(3,:)=outlier_remover3(tt(3,1:S:M),thresh);
%                 tt_out(3,:)=outlier_remover2(tt(3,1:S:M),K,thresh);
                hold on
                plot((0:S:M-1)*dt,tt_out(3,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,1:S:M,3));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S1R3: ','Filtration=',filter,' Resampling=',resam];
            U.S=1;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 4 % S1R4
            figure;
            plot((0:S:M-1)*dt,tt(4,1:S:M),'o');
            if outlierFlg
%                 tt_out(4,:)=outlier_remover(tt(4,1:S:M),0.15);
                tt_out(4,:)=outlier_remover3(tt(4,1:S:M),thresh);
%                 tt_out(4,:)=outlier_remover2(tt(4,1:S:M),K,thresh);
                hold on
                plot((0:S:M-1)*dt,tt_out(4,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,1:S:M,4));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S1R4: ','Filtration=',filter,' Resampling=',resam];
            U.S=1;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 5 % S1R5
            figure;
            plot((0:S:M-1)*dt,tt(5,1:S:M),'o');
            if outlierFlg
%                 tt_out(5,:)=outlier_remover(tt(5,1:S:M),0.15);
                tt_out(5,:)=outlier_remover3(tt(5,1:S:M),thresh);
%                 tt_out(5,:)=outlier_remover2(tt(5,1:S:M),K,thresh);
                hold on
                plot((0:S:M-1)*dt,tt_out(5,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,1:S:M,5));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S1R5: ','Filtration=',filter,' Resampling=',resam];
            U.S=1;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 6 % S2R1
            figure;
            plot((1:S:M)*dt,tt(6,2:S:M+1),'o');
            if outlierFlg
%                 tt_out(6,:)=outlier_remover(tt(6,2:S:M+1),0.15);
                tt_out(6,:)=outlier_remover3(tt(6,2:S:M+1),thresh);
%                 tt_out(6,:)=outlier_remover2(tt(6,2:S:M+1),K,thresh);
                hold on
                plot((1:S:M)*dt,tt_out(6,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,2:S:M+1,1));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S2R1: ','Filtration=',filter,' Resampling=',resam];
            U.S=2;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 7 % S2R2
            figure;
            plot((1:S:M)*dt,tt(7,2:S:M+1),'o');
            if outlierFlg
%                 tt_out(7,:)=outlier_remover(tt(7,2:S:M+1),0.15);
                tt_out(7,:)=outlier_remover3(tt(7,2:S:M+1),thresh);
%                 tt_out(7,:)=outlier_remover2(tt(7,2:S:M+1),K,thresh);
                hold on
                plot((1:S:M)*dt,tt_out(7,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,2:S:M+1,2));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S2R2: ','Filtration=',filter,' Resampling=',resam];
            U.S=2;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 8 % S2R3
            figure;
            plot((1:S:M)*dt,tt(8,2:S:M+1),'o');
            if outlierFlg
%                 tt_out(8,:)=outlier_remover(tt(8,2:S:M+1),0.15);
                tt_out(8,:)=outlier_remover3(tt(8,2:S:M+1),thresh);
%                 tt_out(8,:)=outlier_remover2(tt(8,2:S:M+1),K,thresh);
                hold on
                plot((1:S:M)*dt,tt_out(8,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,2:S:M+1,3));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S2R3: ','Filtration=',filter,' Resampling=',resam];
            U.S=2;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 9 % S2R4
            figure;
            plot((1:S:M)*dt,tt(9,2:S:M+1),'o');
            if outlierFlg
%                 tt_out(9,:)=outlier_remover(tt(9,2:S:M+1),0.15);
                tt_out(9,:)=outlier_remover3(tt(9,2:S:M+1),thresh);
%                 tt_out(9,:)=outlier_remover2(tt(9,2:S:M+1),K,thresh);
                hold on
                plot((1:S:M)*dt,tt_out(9,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,2:S:M+1,4));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S2R4: ','Filtration=',filter,' Resampling=',resam];
            U.S=2;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 10 % S2R5
            figure;
            plot((1:S:M)*dt,tt(10,2:S:M+1),'o');
            if outlierFlg
%                 tt_out(10,:)=outlier_remover(tt(10,2:S:M+1),0.15);
                tt_out(10,:)=outlier_remover3(tt(10,2:S:M+1),thresh);
%                 tt_out(10,:)=outlier_remover2(tt(10,2:S:M+1),K,thresh);
                hold on
                plot((1:S:M)*dt,tt_out(10,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,2:S:M+1,5));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S2R5: ','Filtration=',filter,' Resampling=',resam];
            U.S=2;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 11 % S3R1
            figure;
            plot((2:S:M+1)*dt,tt(11,3:S:M+2),'o');
            if outlierFlg
%                 tt_out(11,:)=outlier_remover(tt(11,3:S:M+2),0.15);
                tt_out(11,:)=outlier_remover3(tt(11,3:S:M+2),thresh);
%                 tt_out(11,:)=outlier_remover2(tt(11,3:S:M+2),K,thresh);
                hold on
                plot((2:S:M+1)*dt,tt_out(11,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,3:S:M+2,1));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S3R1: ','Filtration=',filter,' Resampling=',resam];
            U.S=3;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 12 % S3R2
            figure;
            plot((2:S:M+1)*dt,tt(12,3:S:M+2),'o');
            if outlierFlg
%                 tt_out(12,:)=outlier_remover(tt(12,3:S:M+2),0.15);
                tt_out(12,:)=outlier_remover3(tt(12,3:S:M+2),thresh);
%                 tt_out(12,:)=outlier_remover2(tt(12,3:S:M+2),K,thresh);
                hold on
                plot((2:S:M+1)*dt,tt_out(12,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,3:S:M+2,2));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S3R2: ','Filtration=',filter,' Resampling=',resam];
            U.S=3;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 13 % S3R3
            figure;
            plot((2:S:M+1)*dt,tt(13,3:S:M+2),'o');
            if outlierFlg
%                 tt_out(13,:)=outlier_remover(tt(13,3:S:M+2),0.15);
                tt_out(13,:)=outlier_remover3(tt(13,3:S:M+2),thresh);
%                 tt_out(13,:)=outlier_remover2(tt(13,3:S:M+2),K,thresh);
                hold on
                plot((2:S:M+1)*dt,tt_out(13,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,3:S:M+2,3));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S3R3: ','Filtration=',filter,' Resampling=',resam];
            U.S=3;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 14 % S3R4
            figure;
            plot((2:S:M+1)*dt,tt(14,3:S:M+2),'o');
            if outlierFlg
%                 tt_out(14,:)=outlier_remover(tt(14,3:S:M+2),0.15);
                tt_out(14,:)=outlier_remover3(tt(14,3:S:M+2),thresh);
%                 tt_out(14,:)=outlier_remover2(tt(14,3:S:M+2),K,thresh);
                hold on
                plot((2:S:M+1)*dt,tt_out(14,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,3:S:M+2,4));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S3R4: ','Filtration=',filter,' Resampling=',resam];
            U.S=3;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 15 % S3R5
            figure;
            plot((2:S:M+1)*dt,tt(15,3:S:M+2),'o');
            if outlierFlg
%                 tt_out(15,:)=outlier_remover(tt(15,3:S:M+2),0.15);
                tt_out(15,:)=outlier_remover3(tt(15,3:S:M+2),thresh);
%                 tt_out(15,:)=outlier_remover2(tt(15,3:S:M+2),K,thresh);
                hold on
                plot((2:S:M+1)*dt,tt_out(15,:),'*r');
            end
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            U.s=squeeze(s_r(:,3:S:M+2,5));
            U.dts=dts;
            U.dt=dt*S;
            U.title=['S3R5: ','Filtration=',filter,' Resampling=',resam];
            U.S=3;
            set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
    end
end



