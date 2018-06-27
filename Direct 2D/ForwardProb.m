function [tt_pert, dist, azi] = ForwardProb(xsrc, ysrc, xrcv, yrcv, xaxis, yaxis, c, u, v, c0, sigt)
%FORWARDPROB Solves the forward problem for travel-time tomography.
%
%Format:
%  [tt_pert, dist, azi] = ForwardProb(xsrc, ysrc, xrcv, yrcv, xaxis, yaxis, c, u, v, c0, sigt);
%
%Input:
%(xsrc,ysrc) = Ns x 1 vectors of the x,y coordinates of the source positions, where Ns is the number of sources
%(xrcv,yrcv) = Nr x 1 vectors of the x,y coordinates of the rcvr positions, where Nr is the number of rcvrs
%(xaxis,yaxis) = x,y axes along which the atmospheric fields are specified
%(c,u,v) = atmospheric sound speed and wind components on 2-D grid specified by (xaxis,yaxis) 
%         (u,v are optional, default to zero)
%c0 = reference sound speed (optional, default to mean of c)
%sigt = a struct variable sigt.n containing the standard deviation of travel-time from measurement noise (optional,
%defaults to zero); 
%
%Output:
%tt_pert = Ns x Nr array of travel time perturbations
%dist = lengths of ray paths
%azi = angles of ray paths

% Make sure source/receiver coordinates are column vectors.
xsrc = xsrc(:);   
ysrc = ysrc(:);   
xrcv = xrcv(:);   
yrcv = yrcv(:);
Ns = length(xsrc);
Nr = length(xrcv);

% Set default arguments if not present.
if ~exist('u'),
  u = zeros(size(c));
end
if ~exist('v'),
  v = zeros(size(c));
end
if ~exist('c0'),
  c0 = mean(c(:));
end
if ~exist('sigt'),
  sig_t = 0;
else
    sig_t=sigt.n;
end

% Determine path distances and azimuths.
azi = atan2((ones(Ns,1)*yrcv'-ysrc*ones(1,Nr)), ...
  (ones(Ns,1)*xrcv'-xsrc*ones(1,Nr)));
dist = sqrt((ones(Ns,1)*yrcv'-ysrc*ones(1,Nr)).^2 + ...
  ((ones(Ns,1)*xrcv'-xsrc*ones(1,Nr))).^2);

% Loop through the propagation paths, determining travel time on each one.
tt_pert = zeros(Ns, Nr);
for m = 1:Ns,
  for n = 1:Nr,
    
    % Determine the propagation angle, distance, and effective sound-speed perturbation along the path.
    integ = -(c-c0 + u*cos(azi(m,n)) + v*sin(azi(m,n)))/c0^2;
    
    % Determine the travel time along the path.
    if dist(m,n)==0
        tt_pert(m,n)=0;
    else
        tt_pert(m,n) = quad0(@path_interp, 0, dist(m,n), 1e-8, 0, ...
            xsrc(m), ysrc(m), azi(m,n), integ, xaxis, yaxis);
    end
  end
end

% Include noise contribution.
tt_pert = tt_pert + sig_t*randn(Ns,Nr);
