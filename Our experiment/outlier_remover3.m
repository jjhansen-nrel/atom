function tt=outlier_remover3(t,thresh)
% Input:
% t - a matrix of original travel times, e. g., [15,40];
% tresh - a treshold for large time differences, e.g., tresh=0.15 means
% that 15% of the largest time differences will be suspected as possible
% outliers. The bigger this parameter, the smoother will be filtered travel
% times.
% Output:
% tt - a matrix of filtered travel times of the same size as t.

[Nrays,Ntime]=size(t);
tt=zeros(Nrays,Ntime);
% determine cutt-off limit for Q1
j1=fix(0.25*Ntime);
% determine cutt-off limit for Q3
j3=fix(0.75*Ntime);
% for the median
j2=fix(0.5*Ntime);
outlierFlg=1;
outliers0=[];
K=8;
for i=1:Nrays
%     while outlierFlg 
        outlierFlg=0;
        a=t(i,:);
        % determine the residuals
%         [s,y]=minabs(a);
%         [s,y]=minFFTres(a);
        y=mov_av(a,K,[],'full');
        y=y.';
        s=a-y;
        % find Q3
        s1=sort(s);
        Q3=s1(j3);
        Q2=s1(j2);
        Q1=s1(j1);
        iqb=Q3-Q1;
        % outlier indices
        outliers=find(s > Q2+thresh*iqb | s < Q2-thresh*iqb);
        if ~isempty(outliers) && ~isequal(outliers0,outliers)
            validInd=setdiff(1:Ntime,outliers);
            a(outliers)=interp1(validInd,a(validInd),outliers,'linear','extrap');
            outlierFlg=1;
            outliers0=outliers;
        end
        t(i,:)=a;
%     end
    tt(i,:)=a;
end

function [minRes,y]=minabs(x)
lx=length(x);
x0=1:lx;
lae=inf;
for i=1:lx-1
    for j=i+1:lx
        a=(x(j)-x(i))/(j-i);
        y1=x(i)+a*(x0-i);
        res=x-y1;
        e=sum(abs(res));
        if e < lae
            lae=e;
            minRes=res;
            y=y1;
        end
    end
end

function [s,y]=minFFTres(a)
% ma=median(a);
ma=mean(a);
la=length(a);
la2=round(la/2);
a1=[zeros(1,la2),a-ma,zeros(1,la-la2)];
la1=length(a1);
 fs=fft(a1);
 i0=round(10*la1/la);
 i1=round(30*la1/la);
 fs(i0:i1)=0;
 a2=real(ifft(fs));
 y=a2(la2+1:la2+la)+ma;
 s=a-y;
 
