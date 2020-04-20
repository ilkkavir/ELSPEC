function [h,ts,te,pp,ppstd,loc] = readGUISDAPpp( ppdir , exp , radar , ...
                                             version , tres )
%
% Read GUISDAP raw electron densities (power profiles) and their
% standard deviations from files.
%
% To calculate sufficient data with GUISDAP, use the options
% analysis_ppshortlags=1
% analysis_pponly=1
%
% INPUT:
%  ppdir   path to the data directory
%  exp     name of the eiscat experiment
%  radar   name of the radar 'uhf', 'vhf' or 'esr'
%  version experiment version number
%  tres    "type" of time resolution 'best' or 'dump'
%
% Currently available combinations of exp,radar,version are
%  'beata','u',1
%  'arc1','u',1
%  'arc1','u',2
%
% OUTPUT:
%  h     heights [km]
%  ts    integration start times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  te    integration end times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  pp    a matrix of electron density values [m^-3]
%  ppstd standard deviations of the electron densities [m^-3]
%  loc     [latitude (deg north) longitude (deg east) height (km)]
%          of the radar site
%
% IV 2016
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% list mat files
fp = listGUISDAPfiles( ppdir );

% if no files available
if isempty(fp)
    h=[];
    ts=[];
    te=[];
    pp=[];
    ppstd=[];
    loc=[];
    return
end

switch lower(exp)
    
  case 'beata'
    switch lower(radar)
      case 'uhf'
        [h,ts,te,pp,ppstd,loc,azel] = readGUISDAPpp_beata_uhf( fp , version ...
                                                      , tres );
      otherwise
        error(['Experiment ' exp ' is not implemented for radar ' radar])
    end
  case 'arc1'
    switch lower(radar)
      case 'uhf'
        [h,ts,te,pp,ppstd,loc,azel] = readGUISDAPpp_arc1_uhf( fp , version ...
                                                      , tres );
      otherwise
        error(['Experiment ' exp ' is not implemented for radar ' radar])
    end
  otherwise
    error(['Experiment ' exp ' is not implemented in readGUISDAPpp']);
end

end



function [h,ts,te,pp,ppstd,loc,azel] = readGUISDAPpp_beata_uhf( ff , version ...
                                                      , tres )
%
% [h,ts,te,pp,ppstd,loc,azel] = readGUISDAPpp_beata_uhf( ff , version )
%
% Read raw electron densities from GUISDAP UHF beata result files.
% Sufficient power profiles are available ONLY IF the analysis was
% run with parameters analysis_ppshortlags=1 and analysis_pponly=1
%
% This function is based on guisdap_pp_reader by Björn Gustavsson.
%
% INPUT:
%   ff       a char vector of file names
%   version  exp version number, ignored at the moment
%   tres     "type" of time resolution 'best' or 'dump'
%
% OUTPUT:
%  h     heights [km]
%  ts    integration start times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  te    integration end times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  pp    a matrix of electron density values [m^-3]
%  ppstd standard deviations of the electron densities [m^-3]
%  loc   [latitude (deg north) longitude (deg east) height (km)]
%        of the radar site
%  azel  azimuth and elevation angles of the radar beam
%
% IV 2016

% load the first file to get the matrix sizes
load(ff{1});
i_end = find(diff(r_pprange)<0);
nh = i_end(1);
nf = length(ff);
loc = r_XMITloc;

ppmat = NaN(nh,length(i_end)+1);

% note to self: this is the only place where both options are
% checked. After this if-statement, we asssume that tres is 'best'
% if it is not 'dump'. Remember to change the other if-statements
% if other options are implemented.
if lower(tres)=='dump'
    pp = NaN(nh,nf);
    ppstd = NaN(nh,nf);
    ts = NaN(nf,1);
    te = NaN(nf,1);
    h = NaN(nh,nf);
    azel = NaN(nf,2);
elseif lower(tres)=='best'
    nppdump = (length(i_end)+1);
    npp = nf * nppdump;
    pp = NaN(nh,npp);
    ppstd = NaN(nh,npp);
    ts = NaN(npp,1);
    te = NaN(npp,1);
    h = NaN(nh,npp);
    azel = NaN(npp,2);
else
    error(['Unknown time resolution type',tres])
end


for k = 1:length(ff)

    load(ff{k})

    % there are several profiles in the file, find their ends
    i_end = find(diff(r_pprange)<0);

    % start indices of profiles
    i2 = [1;i_end+1];

    % end indices of profiles
    i3 = [i_end;length(r_pp)];

    % start and end times of integration, converted into unix time
    if lower(tres) == 'dump'
        ts(k) = date2unixtime(r_time(1,:));
        te(k) = date2unixtime(r_time(2,:));
    else
        ts( (k-1)*nppdump+1 ) = date2unixtime(r_time(1,:));
        te( k*nppdump ) = date2unixtime(r_time(2,:));
        dtpp = (te( k*nppdump ) - ts( (k-1)*nppdump+1 ) ) / nppdump;
        for kpp = 2:nppdump
            ts( (k-1)*nppdump+kpp ) = ts( (k-1)*nppdump+1 ) + (kpp-1)*dtpp;
            te( k*nppdump-kpp+1 ) = te( k*nppdump ) - (kpp-1)*dtpp;
        end
    end

    % heights and pointing directions
    if lower(tres) == 'dump'
        h( 1 : (i3(1)-i2(1)+1) , k)  = r_pprange(i2(1):i3(1))*sin(r_el*pi/180);
        azel(k,:) = [r_az r_el];
    else
        h( 1 : (i3(1)-i2(1)+1 ) , (k-1)*nppdump+1 ) = r_pprange(i2(1):i3(1))*sin(r_el*pi/ ...
                                                          180);
        azel((k-1)*nppdump+1,:) = [r_az r_el];
        for kpp=2:nppdump
            h( : , (k-1)*nppdump+kpp ) = h( : , (k-1)*nppdump+1 );
            azel((k-1)*nppdump+kpp) = [r_az r_el];
        end
    end

    % read the power profiles
    for i4 = 1:length(i2),
        ppmat( 1 : (i3(i4)-i2(i4)+1) ,i4) = r_pp(i2(i4):i3(i4));
    end
    ppsize = size(ppmat,1);

    if lower(tres)=='dump'
        % integrate the profiles together
        pp(1:ppsize,k) = mean(ppmat,2);
        ppstd(1:ppsize,k) = std(ppmat,0,2)./sqrt(size(ppmat,2));
    else
        % keep the profiles separate
        pp( 1 : ppsize , ((k-1)*nppdump+1) : (k*nppdump) ) = ppmat;
        ppstd( 1 : ppsize , (k-1)*nppdump+1 ) = std(ppmat,0,2);
        for kpp = 2:nppdump
            ppstd( 1 : ppsize , (k-1)*nppdump+kpp ) = ppstd( : , (k-1)*nppdump+1 );
        end
    end

end

% if MATLAB changed the matrix sizes, there are some zeros that
% must be replaced with NaN's
naninds = h==0;
h(naninds) = NaN;
ppstd(naninds) = NaN;
pp(naninds) = NaN;


% smoother variances...
if lower(tres)=='best'
    for ppind=1:npp
        ppstd(:,ppind) = std(pp(:,max(1,ppind-5):min(npp,ppind+5))');
    end
end

end





function [h,ts,te,pp,ppstd,loc,azel] = readGUISDAPpp_arc1_uhf( ff , version ...
                                                      , tres )
%
% [h,ts,te,pp,ppstd,loc] = readGUISDAPpp_arc1_uhf( ff , version )
%
% Read raw electron densities from GUISDAP UHF arc1 result files.
% Sufficient power profiles are available ONLY IF the analysis was
% run with parameters analysis_ppshortlags=1 and analysis_pponly=1
%
% This function is based on guisdap_pp_reader by Björn Gustavsson.
%
% INPUT:
%   ff       a char vector of file names
%   version  exp version number, ignored at the moment
%   tres     "type" of time resolution 'best' or 'dump'
%
% OUTPUT:
%  h     heights [km]
%  ts    integration start times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  te    integration end times, unix time (seconds since 1970-01-01 00:00:00 UT)
%  pp    a matrix of electron density values [m^-3]
%  ppstd standard deviations of the electron densities [m^-3]
%  loc   [latitude (deg north) longitude (deg east) height (km)]
%        of the radar site
%  azel  beam azimuth and elevation angles
%
% IV 2016

switch version
  case 1
    nprofstep = 3;
  case 2
    nprofstep = 2;
  otherwise
    error('unknown version number')
end

% load the first file to get the matrix sizes
load(ff{1});
i_end = find(diff(r_pprange)<0);
i_end = i_end(1:nprofstep:end); % there are four profiles, we want only the
                        % first one...
nh = i_end(1);
nf = length(ff);

loc = r_XMITloc;

ppmat = NaN(nh,length(i_end));

% note to self: this is the only place where both options are
% checked. After this if-statement, we asssume that tres is 'best'
% if it is not 'dump'. Remember to change the other if-statements
% if other options are implemented.
if lower(tres)=='dump'
    pp = NaN(nh,nf);
    ppstd = NaN(nh,nf);
    ts = NaN(nf,1);
    te = NaN(nf,1);
    h = NaN(nh,nf);
    azel = NaN(nf,2);
elseif lower(tres)=='best'
    nppdump = (length(i_end));
    npp = nf * nppdump;
    pp = NaN(nh,npp);
    ppstd = NaN(nh,npp);
    ts = NaN(npp,1);
    te = NaN(npp,1);
    h = NaN(nh,npp);
    azel = NaN(npp,2);
else
    error(['Unknown time resolution type',tres])
end


for k = 1:length(ff)

    load(ff{k})

    % there are several profiles in the file, find their ends
    i_end = find(diff(r_pprange)<0);

    % start indices of profiles
    i2 = [1;i_end(nprofstep:nprofstep:end)+1];

    % end indices of profiles
    i3 = [i_end(1:nprofstep:end)];%length(r_pp)];

    % start and end times of integration, converted into unix time
    if lower(tres) == 'dump'
        ts(k) = date2unixtime(r_time(1,:));
        te(k) = date2unixtime(r_time(2,:));
    else
        ts( (k-1)*nppdump+1 ) = date2unixtime(r_time(1,:));
        te( k*nppdump ) = date2unixtime(r_time(2,:));
        dtpp = (te( k*nppdump ) - ts( (k-1)*nppdump+1 ) ) / nppdump;
        for kpp = 2:nppdump
            ts( (k-1)*nppdump+kpp ) = ts( (k-1)*nppdump+1 ) + (kpp-1)*dtpp;
            te( k*nppdump-kpp+1 ) = te( k*nppdump ) - (kpp-1)*dtpp;
        end
    end

    % heights and beam directions
    if lower(tres) == 'dump'
        h(:,k)  = r_pprange(i2(1):i3(1))*sin(r_el*pi/180);
        azel(k,:) = [r_az r_el];
    else
        h( : , (k-1)*nppdump+1 ) = r_pprange(i2(1):i3(1))*sin(r_el*pi/ ...
                                                          180);
        azel((k-1)*nppdump+1) = [r_az r_el];
        for kpp=2:nppdump
            h( : , (k-1)*nppdump+kpp ) = h( : , (k-1)*nppdump+1 );
            azel((k-1)*nppdump+kpp) = [r_az r_el];
        end
    end

    % read the power profiles
    for i4 = 1:length(i2),
        ppmat(:,i4) = r_pp(i2(i4):i3(i4));
    end

    if lower(tres)=='dump'
        % integrate the profiles together
        pp(:,k) = mean(ppmat,2);
        ppstd(:,k) = std(ppmat,0,2)./sqrt(size(ppmat,2));
    else
        % keep the profiles separate
        pp( : , ((k-1)*nppdump+1) : (k*nppdump) ) = ppmat;
        ppstd( : , (k-1)*nppdump+1 ) = std(ppmat,0,2);
        for kpp = 2:nppdump
            ppstd( : , (k-1)*nppdump+kpp ) = ppstd( : , (k-1)*nppdump+1 );
        end
    end

end

% smoother variances...
if lower(tres)=='best'
    for ppind=1:npp
        ppstd(:,ppind) = std(pp(:,max(1,ppind-5):min(npp,ppind+5))');
    end
end

end