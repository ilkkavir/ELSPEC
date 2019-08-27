function [h,ts,te,par,parstd,loc] = readGUISDAPpar( fitdir )
%
% Read GUISDAP plasma parameter fit results from mat-files
%
% [h,ts,te,par,parstd] = readGUISDAPpar( fitdir )
%
% INPUT:
%  fitdir  path to GUISDAP fit results
%
% OUTPUT:
%  h       heights [km]
%  ts      integration period start times [unixtime]
%  te      integration period end times [unixtime]
%  par     lenth(h)x4xlength(ts) array of plasma parameters
%  parstd  standard deviations of the plasma parameters
%  loc     [latitude (deg north) , longitude (deg east) , altitude
%           (km)] of the radar site
%
% The four plasma parameters are Ne [m^-3], Ti [K], Te [k], and Vi
% [ms^-1]. Failed iterations and results with chi-squared larger
% than 10 are replaced with NaNs.
%
% IV 2016 - 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% list GUISDAP output files
fp = listGUISDAPfiles( fitdir );

% if no files available
if isempty(fp)
    h=[];
    ts=[];
    te=[];
    par=[];
    parstd=[];
    loc=[];
    return
end

% number of files
nf = length(fp);

% load the first file to get initial matrix dimensions
load(fp{1});
nh = length(r_h);

loc = r_XMITloc;

h = NaN(nh,nf);
ts = NaN(nf,1);
te = NaN(nf,1);
par = NaN(nh,4,nf);
parstd = NaN(nh,4,nf);

for k=1:nf

    load(fp{k});

    % start and end times of integration, converted into unix time
    ts(k) = date2unixtime(r_time(1,:));
    te(k) = date2unixtime(r_time(2,:));

    % sometimes we do not have fit from all gates
    nh1 = length(r_h);
    if nh1 >= nh
        nh = nh1;
        ih = 1:nh;
    else
        ih=[];
        for ii = 1:length(r_h)
            [dummy,ih(ii)] = min(abs(h(:,1)-r_h(ii)));
        end
    end

    h(ih,k) = r_h;
    par(ih,:,k) = r_param(:,[1 2 3 5]);
    parstd(ih,:,k) = r_error(:,[1 2 3 5]);

    % remove failed iterations
    %    failed = (r_status ~= 0) | ( r_res(:,1) > 10 );
    %    failed = (r_status ~= 0) | ( r_res(:,1) > 30 );
    failed = (all(r_status ~= [0 3])) | ( r_res(:,1) > 3000 );
    par(failed,:,k) = NaN;
    parstd(failed,:,k) = NaN;

end

par(:,3,:) = par(:,3,:) .* par(:,2,:);
parstd(:,3,:) = parstd(:,2,:).*par(:,3,:) + parstd(:,3,:).*par(:,2,:);