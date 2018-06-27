function path=ttFiltr_r(S,R,xy,t,Lo,Hi,showFig)
% ttFiltr_r is a plausible filtr for travel times for reciprocal
% transmission arrays of transducers; it does the same job as ttFiltr for
% non reciprocal arrays, but hte output format is different.
% Syntax:
% path=ttFiltr(S,R,xy,t,Lo,Hi,showFig);
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
% path - a struct vector of the valid travel paths of the length less or equal
%   to S*R (can be less if some of the travel times have been filtered out) with fields:
%       length - the length of a valid travel path;
%       t - its travel time (sec);
%       index - its index;
%       theta - sx and sy coordinates of the unit vector s in the direction of
%           the ray propagation;
%       xy0 - the x and y coordinates of the valid ray's source;
%       xy - the x and y coordinates of the valid ray's receiver.

if length(t)~=S*R
    error('Length of travel time vector is not equal to number of rays.');
end

for i=1:S % loop over all sources
    for j=1:R % loop over all receivers
        ind=R*(i-1)+j;
        b=xy(S+j,1)-xy(i,1); % difference of x between the j-th receiver and i-th source 
        a=xy(S+j,2)-xy(i,2); % difference of y between the j-th receiver and i-th source
        lLi=sqrt(a*a+b*b);
        if lLi==0
            path(ind).theta=[NaN NaN];
        else
            Sx=b/lLi; % x coordinate of the unit vector S - direction of the group velocity u
            Sy=a/lLi; % y coordinate of the unit vector S
            path(ind).theta=[Sx Sy];
        end
        path(ind).xy0=[xy(i,1) xy(i,2)];
        path(ind).xy=[xy(S+j,1) xy(S+j,2)];
        path(ind).index=ind;
        path(ind).length=lLi;
        path(ind).t=t(ind);
    end
end
i=isnan(t) | t==0 | [path.length]'==0;
path(i)=[];
c=[path.length]'./[path.t]';
if showFig
    figure;
    plot([path.index]',c,'*');
    hold on
    plot([path(1).index path(end).index],[Lo Lo],'r','linewidth',2);
    plot([path(1).index path(end).index],[Hi Hi],'r','linewidth',2);
    xlabel('Index of ray','fontSize',12,'fontweight','bold');
    ylabel('Aver. speed along ray (m/s)','fontSize',12,'fontweight','bold');
end
i= c<Lo | c>Hi;
path(i)=[];
