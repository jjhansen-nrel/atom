function [Li,t,in,theta]=ttFiltr(S,R,xy,t,Lo,Hi,showFig)
% ttFiltr is a plausible filtr for travel times;
% Syntax:
% [Li,t,in,theta]=ttFiltr(S,R,xy,t,Lo,Hi,showFig);
% Inputs:
% S - number of sound sources, integer;
% R - number of receivers, integer;
% xy - the x and y coordinated of sources and receivers; a
%   matrix [S+R, 2]; the first S rows contain the coordinates of sources,
%   the last R rows contain the coordinates of receivers;
% t - travel times (sec), a vector of length S*R;
% Lo - a lower limit for the filtration; travel times which corresponds to
%   the effective sound speed lower then Lo will be ommited;
% Hi - an upper limit for the filtration; travel times which corresponds to
%   the effective sound speed higher then Hi will be ommited;
% showFig - a flag to display a graph; if showFig=1, the graph showing Lo,
%   Hi and effective sound speed will be displayed, otherwise, it will not be
%   displayed.
% Outputs:
% Li - the lengths of the valid travel paths, a vector of the length less or equal
%   to S*R (can be less if some of the travel times are filtered out);
% t - filtered travel times (sec), a vector of the same size as Li;
% in - indeces of valid travel times, a vector of the same size as Li;
% theta - sx and sy coordinates of the unit vector s in the direction of
%   the ray propagation, a matrix [number of valid rays, 2].

if length(t)~=S*R
    error('Length of tr. times mismatches to the number of paths.');
end
Li=zeros(S*R,1);
theta=zeros(S*R,2);
for i=1:S % loop over all sources
    for j=1:R % loop over all receivers
        b=xy(S+j,1)-xy(i,1); % difference of x between the j-th receiver and i-th source 
        a=xy(S+j,2)-xy(i,2); % difference of y between the j-th receiver and i-th source
        lLi=sqrt(a*a+b*b);
        Li(R*(i-1)+j)=lLi;
        if lLi~=0
            Sx=b/lLi; % x coordinate of the unit vector S - direction of the group velocity u
            Sy=a/lLi; % y coordinate of the unit vector S
        else
            Sx=NaN; % x coordinate of the unit vector S - direction of the group velocity u
            Sy=NaN; % y coordinate of the unit vector S
        end
            theta(R*(i-1)+j,:)=[Sx Sy];
        
    end
end
in=[1:length(t)]';
i=isnan(t) | t==0 | Li==0 | t<5e-3;
t(i)=[];
Li(i)=[];
in(i)=[];
theta(i,:)=[];
c=Li./t;
if showFig
    figure;
    plot(in,c,'*');
    hold on
    plot([in(1) in(end)],[Lo Lo],'r','linewidth',2);
    plot([in(1) in(end)],[Hi Hi],'r','linewidth',2);
    xlabel('Path index','fontSize',12,'fontweight','bold');
    ylabel('Aver. speed along the path (m/s)','fontSize',12,'fontweight','bold');
end
i=find(c<Lo | c>Hi);
t(i)=[];
Li(i)=[];
in(i)=[];
theta(i,:)=[];