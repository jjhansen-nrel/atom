function outlier_remover_m_batch(use_t_out_flg)
if nargin==0
    use_t_out_flg=1; % load most recent travel times, already filtered before, and start from here.
end
path_name=uigetdir('Pick a folder with data');

if path_name == 0
    return
end
path_name=[path_name,'\'];
F=dir([path_name,'*.mat']);
% dt=60/M; % record interval in sec.
outlierFlg=1;
tt0=[];
R=5;
file_name=F(1).name;
if isequal('tt_out',file_name(1:6))
    load([path_name,file_name]);
    [Nrays,Ntime]=size(tt_out); %#ok<*NODEF>
    S=3;
    M=Ntime*3;
    tt=zeros(Nrays,M);
    resam='1';
    filter='1';
    filename=file_name(8:end-4);
    tt(1:R,1:S:end)=tt_out(1:R,:); %#ok<NODEF>
    tt(R+1:2*R,2:S:end)=tt_out(R+1:2*R,:);
    tt(2*R+1:3*R,3:S:end)=tt_out(2*R+1:3*R,:);
%     dt=60/M; %#ok<NODEF>
else
    for fi=1:length(F)
        file_name=F(fi).name;
        load([path_name,file_name]);
        filename=[file_name(4:7),'/',file_name(8:9),'/',file_name(10:11),' ',file_name(12:13),':',file_name(14:15),':',file_name(16:17)];
        if fi==1
            filename1=file_name;
        elseif fi==length(F)
            filename2=file_name;
        end
        resam=num2str(upSampleFlg);
        filter=num2str(filterFlg);
        
        if exist('tt_out','var')
            if use_t_out_flg
                tt(1:R,1:S:end)=tt_out(1:R,:); %#ok<NODEF>
                tt(R+1:2*R,2:S:end)=tt_out(R+1:2*R,:);
                tt(2*R+1:3*R,3:S:end)=tt_out(2*R+1:3*R,:);
            end
        end
        tt0=[tt0 tt];
    end
%     dt=60/M; %#ok<NODEF>
    M=M*length(F);
    Ntime=M/S;
    tt_out=zeros(size(tt,1),Ntime);
    tt=tt0;
end
Choice=1;
dt=0.5;

while Choice~=16
    Choice=menu(filename,'S1R1','S1R2','S1R3','S1R4', 'S1R5',...
        'S2R1','S2R2','S2R3','S2R4', 'S2R5',...
        'S3R1','S3R2','S3R3','S3R4', 'S3R5','Exit');
    switch Choice
        case 1 % S1R1
            figure;
            a=tt(1,1:S:M);
            plot((0:S:M-1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            % determine outliers
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round(x/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(1,:)=a;
                hold on
                plot((0:S:M-1)*dt,tt_out(1,:),'*r');
            end
            %             U.s=squeeze(s_r(:,1:S:M,1));
            %             U.dts=dts;
            %             U.dt=dt*S;
            %             U.title=['S1R1: ','Filtration=',filter,' Resampling=',resam];
            %             U.S=1;
            %             set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 2 % S1R2
            a=tt(2,1:S:M);
            figure;
            plot((0:S:M-1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            % determine outliers
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round(x/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(2,:)=a;
                hold on
                plot((0:S:M-1)*dt,tt_out(2,:),'*r');
            end
            
            %                 U.s=squeeze(s_r(:,1:S:M,2));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S1R2: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=1;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 3 % S1R3
            a=tt(3,1:S:M);
            figure;
            plot((0:S:M-1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round(x/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(3,:)=a;
                hold on
                plot((0:S:M-1)*dt,tt_out(3,:),'*r');
            end
            
            %                 U.s=squeeze(s_r(:,1:S:M,3));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S1R3: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=1;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 4 % S1R4
            a=tt(4,1:S:M);
            figure;
            plot((0:S:M-1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round(x/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(4,:)=a;
                hold on
                plot((0:S:M-1)*dt,tt_out(4,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,1:S:M,4));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S1R4: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=1;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 5 % S1R5
            a=tt(5,1:S:M);
            figure;
            plot((0:S:M-1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S1R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round(x/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(5,:)=a;
                hold on
                plot((0:S:M-1)*dt,tt_out(5,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,1:S:M,5));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S1R5: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=1;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 6 % S2R1
            a=tt(6,2:S:M+1);
            figure;
            plot((1:S:M)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(6,:)=a;
                hold on
                plot((1:S:M)*dt,tt_out(6,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,2:S:M+1,1));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S2R1: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=2;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 7 % S2R2
            a=tt(7,2:S:M+1);
            figure;
            plot((1:S:M)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(7,:)=a;
                hold on
                plot((1:S:M)*dt,tt_out(7,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,2:S:M+1,2));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S2R2: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=2;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 8 % S2R3
            a=tt(8,2:S:M+1);
            figure;
            plot((1:S:M)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(8,:)=a;
                hold on
                plot((1:S:M)*dt,tt_out(8,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,2:S:M+1,3));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S2R3: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=2;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 9 % S2R4
            a=tt(9,2:S:M+1);
            figure;
            plot((1:S:M)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(9,:)=a;
                hold on
                plot((1:S:M)*dt,tt_out(9,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,2:S:M+1,4));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S2R4: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=2;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 10 % S2R5
            a=tt(10,2:S:M+1);
            figure;
            plot((1:S:M)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S2R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(10,:)=a;
                hold on
                plot((1:S:M)*dt,tt_out(10,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,2:S:M+1,5));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S2R5: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=2;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 11 % S3R1
            a=tt(11,3:S:M+2);
            figure;
            plot((2:S:M+1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R1: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-2*dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(11,:)=a;
                hold on
                plot((2:S:M+1)*dt,tt_out(11,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,3:S:M+2,1));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S3R1: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=3;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 12 % S3R2
            a=tt(12,3:S:M+2);
            figure;
            plot((2:S:M+1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R2: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-2*dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(12,:)=a;
                hold on
                plot((2:S:M+1)*dt,tt_out(12,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,3:S:M+2,2));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S3R2: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=3;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 13 % S3R3
            a=tt(13,3:S:M+2);
            figure;
            plot((2:S:M+1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R3: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-2*dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(13,:)=a;
                hold on
                plot((2:S:M+1)*dt,tt_out(13,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,3:S:M+2,3));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S3R3: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=3;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 14 % S3R4
            a=tt(14,3:S:M+2);
            figure;
            plot((2:S:M+1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R4: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-2*dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(14,:)=a;
                hold on
                plot((2:S:M+1)*dt,tt_out(14,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,3:S:M+2,4));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S3R4: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=3;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
        case 15 % S3R5
            a=tt(15,3:S:M+2);
            figure;
            plot((2:S:M+1)*dt,a,'-o');
            xlabel('Time (s)','fontweight','bold','fontsize',12)
            ylabel('Travel time (ms)','fontweight','bold','fontsize',12)
            title_txt=['S3R5: ','Filtration=',filter,' Resampling=',resam];
            title(title_txt,'fontweight','bold','fontsize',12)
            
            if outlierFlg
                button=1;
                outliers=[];
                while button==1
                    [x, ~, button] = ginput(1);
                    k=round((x-2*dt)/(S*dt))+1;
                    if k < 1
                        k=1;
                    elseif k > Ntime
                        k=Ntime;
                    end
                    if button==1
                        hold on
                        i=find(outliers==k, 1);
                        if isempty(i)
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','b')
                            outliers=[outliers k];
                        else
                            plot((k-1)*S*dt+2*dt,a(k),'o','MarkerFaceColor','w')
                            outliers(i)=[];
                        end
                    end
                end
                if ~isempty(outliers)
                    validInd=setdiff(1:Ntime,outliers);
                    a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
                end
                tt_out(15,:)=a;
                hold on
                plot((2:S:M+1)*dt,tt_out(15,:),'*r');
            end
            %                 U.s=squeeze(s_r(:,3:S:M+2,5));
            %                 U.dts=dts;
            %                 U.dt=dt*S;
            %                 U.title=['S3R5: ','Filtration=',filter,' Resampling=',resam];
            %                 U.S=3;
            %                 set(gca,'UserData',U,'ButtonDownFcn','ClickHandler')
    end
end
save('tt_out.mat','tt_out');
