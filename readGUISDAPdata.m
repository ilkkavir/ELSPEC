function [h,ts,te,pp,ppstd,par,parstd,loc] = readGUISDAPdata( ppdir , ...
                                                  fitdir , hmin , ...
                                                  hmax , tmin , tmax ...
                                                  , exp , radar , ...
                                                      version , tres ...
                                                      )
%
% Read GUISDAP raw densities (power profiles) and fit results from
% files in ppdir and fitdir.
%
% [h,ts,te,pp,ppstd,par,parstd] =
%    readGUISDAPdata( ppdir , fitdir , hmin , hmax , tmin , tmax , ...
%    exp , radar , version )
%
%
% INPUT:
%  ppdir   path to GUISDAP power profile data
%  fitdir  path to GUISDAP fit results
%  hmin    minimum height [km]
%  hmax    maximum height [km]
%  tmin    begin time, array [year, month, day, hour , minute ,
%          second], use empty array to read all data
%  tmax    end time, array [year, month, day, hour , minute ,
%          second], use empty array to read all data
%  exp     EISCAT experiment name
%  radar   EISCAT radar name ['U','V']
%  version EISCAT experiment version number [1,2,3,...]
%  tres    "type" of time resolution 'best' or 'dump'
%
%
% OUTPUT:
%  h       heights [km]
%  ts      integration period start times [unixtime]
%  te      integration period end times [unixtime]
%  pp      raw densities, length(h)xlength(ts) matrix [m^-3]
%  ppstd   standard deviations of the raw densities [m^-3]
%  par     lenth(h)x4xlength(ts) array of plasma parameters
%  parstd  standard deviations of the plasma parameters
%
% The four plasma parameters are Ne [m^-3], Ti [K], Te [k], and Vi
% [ms^-1]. Failed iterations and results with chi-squared larger
% than 10 are replaced with NaNs.
%
% IV 2016 - 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% a special case for empty ppdir
if isempty(ppdir)
    [hpar,ts,te,par,parstd,loc] = readGUISDAPpar( fitdir );
    hind = hpar(:,1)>=hmin & hpar(:,1)<=hmax;
    if isempty(tmin)
        t1 = -Inf;
    else
        t1 = date2unixtime(tmin);
    end
    if isempty(tmax)
        t2 = Inf;
    else
        t2 = date2unixtime(tmax);
    end
    it = find((ts >= t1) & (te <= t2));
    ts = ts(it);
    te = te(it);
    h = hpar(hind,1)';
    par = par(hind,:,it);
    parstd = parstd(hind,:,it);
    pp = squeeze(par(:,1,:));
    ppstd = squeeze(parstd(:,1,:));


    % the guisdap error estimates are not good enough for us, try
    % to calculate an ensamble average
    [d1 d2 d3] = size(parstd);
    for i1 = 1:d1
        for i3 = 1:d3
            tmpne = par(max(1,i1-1):min(d1,i1+1),1,max(1,i3-1): ...
                        min(d3,i3+1));
            ppstd(i1,i3) = std(tmpne(:),'omitnan');
        end
    end
    
    return
end

% read power profiles
[hpp,tspp,tepp,pp1,ppstd1,locpp] = readGUISDAPpp( ppdir , exp , radar , ...
                                            version , tres );

% read fit results
[hpar,tspar,tepar,par1,parstd1,locpar] = readGUISDAPpar( fitdir );

% pick the location
loc = locpp;
if isempty(loc)
    loc = locpar;
end
if isempty(loc)
    switch lower(radar)
      case 'uhf'
        loc = [69.583 19.21 0.03];
      case 'vhf'
        loc = [69.583 19.21 0.03];
      case 'esr'
        loc = [78.1530   16.0290    0.4380];
      case '42m'
        loc = [78.1530   16.0290    0.4380];
      case '32m'
        loc = [78.1530   16.0290    0.4380];
      otherwise
        error(['radar location not found from data files and unknown radar site' radar])
    end
end

% now we have data in different grids 

% first pick the most common height step from hpp, this will be our
% height resolution
hstps = diff( hpp );
hstps = hstps(:);
hstps(isnan(hstps)) = [];
hstp = mode( hstps );

% then pick the most common first height, this is where we will
% start from
hbegs = hpp(1,:);
hbegs(isnan(hbegs)) = [];
hbeg = mode( hbegs );

% then the most common largest height
hmaxs = max( hpp );
hmaxs(isnan(hmaxs)) = [];
hend = mode(hmaxs);

% then finally the heights
h = hbeg:hstp:hend;
h = h(h>=hmin);
h = h(h<=hmax);
nh = length(h);

% use the times from power profiles
ts = tspp;
te = tepp;
if isempty(tmin)
    t1 = -Inf;
else
    t1 = date2unixtime(tmin);
end
if isempty(tmax)
    t2 = Inf;
else
    t2 = date2unixtime(tmax);
end
it = find((ts >= t1) & (te <= t2));
ts = ts(it);
te = te(it);
nt = length(te);

% interpolate everything to the height grid h and pick only
% time-slices it
pp = NaN(nh,nt);
ppstd = NaN(nh,nt);
for k=1:nt
    % remove NaN values from pp before interpolation
    ii = ~isnan(hpp(:,it(k)));
    pp(:,k) = interp1(hpp(ii,it(k)),pp1(ii,it(k)),h);
    ppstd(:,k) = interp1(hpp(ii,it(k)),ppstd1(ii,it(k)),h);
end

% the fit results, the fit results can be optionally skipped with
% empty fitdir.
par = NaN(nh,4,nt);
parstd = NaN(nh,4,nt);
if ~isempty(fitdir)
    ntpar = length(tepar);
    par2 = NaN(nh,4,ntpar);
    parstd2 = NaN(nh,4,ntpar);
    for k=1:ntpar
        ii = ~isnan(hpar(:,k));
        % interpolate, some tricks at edges
        par2(:,:,k) = interp1([0;hpar(ii,k);1e6],[par1(ii(1),:,k);par1(ii,:,k);par1(ii(end),:,k)],h,'linear','extrap');
        parstd2(:,:,k) = interp1([0;hpar(ii,k);1e6],[parstd1(ii(1),:,k);parstd1(ii,:,k);parstd1(ii(end),:,k)],h,'linear','extrap');
    end
    for k=1:nh
        for l=1:4
            par(k,l,:) = interp1([0;tepar;1e20],[squeeze(par2(k,l,1));squeeze(par2(k,l,:));squeeze(par2(k,l,end))],te,'linear','extrap');
            parstd(k,l,:) = interp1([0;tepar;1e20],[squeeze(parstd2(k,l,1));squeeze(parstd2(k,l,:));squeeze(parstd2(k,l,end))],te,'linear','extrap');
        end
    end
end



end