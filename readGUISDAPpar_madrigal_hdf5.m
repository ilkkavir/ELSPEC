function [h,ts,te,par,parstd,loc,azel,I] = readGUISDAPpar_madrigal_hdf5(madfile,FAdev)
%
% Read GUISDAP plasma parameter fit results from Madrigal hdf5 files
%
%  [h,ts,te,par,parstd,loc] = readMadrigalPar_hdf5(madfile,FAdev)
%
%  INPUT:
%    madfile   full path to a madrigal hdf5 file
%    FAdev     maximum beam direction deviation from field-aligned [deg]
%
% OUTPUT:
%  h       heights [km]
%  ts      integration period start times [unixtime]
%  te      integration period end times [unixtime]
%  par     lenth(h)x4xlength(ts) array of plasma parameters
%  parstd  standard deviations of the plasma parameters
%  loc     [latitude (deg north) , longitude (deg east) , altitude
%           (km)] of the radar site
%  azel    azimuth and elevation of the radar beam
%  I       magnetic inclination angle (deg)
%
% The four plasma parameters are Ne [m^-3], Ti [K], Te [k], and Vi
% [ms^-1]. Failed iterations and results with chi-squared larger
% than 30 are replaced with NaNs.
%
% IV 2021
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% read the file contents
maddata = h5read(madfile,'/Data/Table Layout');

if ~isfield(maddata,'tr')
   maddata.tr = maddata.te ./ maddata.ti;
end
if ~isfield(maddata,'dtr')
   maddata.dtr = maddata.dte ./ maddata.ti + maddata.te./maddata.ti.^2 .*maddata.dte;
end

% radar location must be checked from the kinst code
% Turns out that this functions works only with EISCAT data, but keeping the codes for future reference
kinst = maddata.kinst(1);
switch kinst
    % Millstone Hill
  case {30,31,32}
    loc = [42.619, 288.51];
    % PFISR
  case 61
    loc = [65.13, 212.529];
    % EISCAT Tromso
  case {72,74}
    loc = [69.583, 19.21];
    % Sonderstrom
  case 80
    loc = [ 67.0, 309.0];
    % RISR
  case {91,92}
    loc = [74.72955, 265.09424];
    % EISCAT Svalbard
  case 95
    loc = [78.09, 16.02];
  otherwise
    error('Unknown kinst')
end



% timestamps as unix time
tsall = maddata.ut1_unix;
teall = maddata.ut2_unix;

% unique start times (each unique time should be one profile)
ts = unique(tsall);
nt = length(ts);
te = ts.*NaN;
for it=1:nt
    % all data points with this start time
    its = tsall==ts(it);

    hh = maddata.gdalt(its);
    if it==1
        nh = length(hh);
        h = NaN(nh,nt);
        par = NaN(nh,4,nt);
        parstd = NaN(nh,4,nt);
        azel = NaN(nt,2);
    end

    % sometimes we do not have fit from all gates
    nh1 = length(hh);
    if nh1 >= nh
        nh = nh1;
        ih = 1:nh;
    else
        ih=[];
        for ii = 1:length(hh)
            [dummy,ih(ii)] = min(abs(h(:,1)-hh(ii)));
        end
    end
    h(ih,it) = hh;

    par(ih,1,it) = maddata.ne(its);
    parstd(ih,1,it) = maddata.dne(its);
    par(ih,2,it) = maddata.ti(its);
    parstd(ih,2,it) = maddata.dti(its);
    par(ih,3,it) = maddata.tr(its);
    parstd(ih,3,it) = maddata.dtr(its);
    if isfield(maddata,'vo')
       par(ih,4,it) = maddata.vo(its);
       parstd(ih,4,it) = maddata.dvo(its);
    end
    

    % remove failed iterations
    res = maddata.chisq(its);
    if isfield(maddata,'gfit')
        status = maddata.gfit(its);
    else
      status = 0;
    end
    failed = ((status ~= 0) & (status ~=3)) | ( res(:,1) > 30 );
    par(failed,:,it) = NaN;
    parstd(failed,:,it) = NaN;

    % azimuth and elevation
    azel(it,:) = [mean(maddata.azm(its)), mean(maddata.elm(its))];
    
    % a common integration end time for all points in the profile
    te(it) = max(teall(its));
end


par(:,3,:) = par(:,3,:) .* par(:,2,:);
parstd(:,3,:) = parstd(:,2,:).*par(:,3,:) + parstd(:,3,:).*par(:,2,:);


% remove other than field-aligned data
%
% local field-aligned direction in E region
[~,~, D, I,~,~,~,~,~,~] = igrfmagm(110000, loc(1), loc(2), decyear(datetime(ts(1),'convertfrom','posixtime')));
FAele = abs(I);
FAaz = D+180;
if I<0
    FAaz = D
end
while FAaz < 0
    FAaz = FAaz + 360;
end
while FAaz > 360
    FAaz = FAaz - 360;
end

% FAdev degree tolerance to allow changes in field-direction...
%rminds = abs(mod(azel(:,1),360) - FAaz) > abs(FAdev) | abs(abs(90-azel(:,2))-(90-FAele)) > abs(FAdev);

FAvec = [cos(FAele*pi/180)*cos(FAaz*pi/180) cos(FAele*pi/180)*sin(FAaz*pi/180) sin(FAele*pi/180)];
BEAMvecs = [cos(azel(:,2)*pi/180).*cos(azel(:,1)*pi/180) cos(azel(:,2)*pi/180).*sin(azel(:,1)*pi/180) sin(azel(:,2)*pi/180)];
rminds = acos(BEAMvecs*FAvec')*180/pi > FAdev;

h(:,rminds) = [];
ts(rminds) = [];
te(rminds) = [];
par(:,:,rminds) = [];
parstd(:,:,rminds) = [];
azel(rminds,:) = [];



% do not accept Ti > 300 K below 100 km altitude or Ti < 50 K anywhere. If Ne < 1e9, keep the electron densities but do not use the temperatures
[dim1 dim2] = size(h);
for i1 = 1:dim1
    for i2 = 1:dim2
        if h(i1,i2) < 100
            if  par(i1,2,i2) > 300
                par(i1,:,i2) = NaN;
                parstd(i1,:,i2) = NaN;
            end
        end
        if par(i1,2,i2) < 50
            par(i1,:,i2) = NaN;
            parstd(i1,:,i2) = NaN;
        end
        if par(i1,1,i2) < 1e9 
            par(i1,2:end,i2) = NaN;
            parstd(i1,2:end,i2) = NaN;
        end
        if par(i1,2,i2) < 0 | par(i1,3,i2) < 0
            par(i1,:,i2) = NaN;
            parstd(i1,:,i2) = NaN;
        end
    end
end

end

