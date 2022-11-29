function [h,ts,te,par,parstd,loc,azel,I] = readGUISDAPpar( fitdir , FAdev )
%
% Read GUISDAP plasma parameter fit results from mat-files
%
% [h,ts,te,par,parstd] = readGUISDAPpar( fitdir , FAdev )
%
% INPUT:
%  fitdir  path to GUISDAP fit results
%  FAdev   maximum beam direction deviation from field-aligned  [deg]
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
% IV 2016 - 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% if this is a madrigal hdf5 file
[~,~,fileextension]=fileparts(fitdir);
if strcmp(fileextension,'.hdf5')
    [h,ts,te,par,parstd,loc,azel,I] = readGUISDAPpar_madrigal_hdf5(fitdir,FAdev);
    return
end

mergedfile = false;
if isdir(fitdir)
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
        azel = [];
        I=[];
        %        return
    end
    
    % number of files
    nf = length(fp);

    if nf>0
        % load the first file to get initial matrix dimensions
        load(fp{1});
    end
else
    mergedfile = true;
    % this must be a merged matlab file
    try
        fnames = sort(fieldnames(matfile(fitdir)));
    catch
        % if cannot list the file contents
        h=[];
        ts=[];
        te=[];
        par=[];
        parstd=[];
        loc=[];
        return
    end

    l = 1;
    flist = [];
    for k=1:length(fnames)
        fname = char(fnames(k));
        if length(fname)==12
            if fname(1:4)=='data'
                flist(l).file = str2num(fname(5:12));
                flist(l).fname = fullfile(fileparts(fitdir),[fname(5:12) '.mat']);
                l = l+1;
            end
        end
    end
    nf = length(flist);
    if nf==0
        % if no files available
        h=[];
        ts=[];
        te=[];
        par=[];
        parstd=[];
        loc=[];
        return
    end
    fname = ['data' flist(1).fname((end-11):(end-4))];
    tmpdata = getfield(load(fitdir,fname),fname);
    filefields = fieldnames(tmpdata);
    for ii=1:length(filefields)
        eval([char(filefields(ii)) '=tmpdata.' char(filefields(ii)) ';']);
    end
end

if nf>0
    
    nh = length(r_h);
    
    loc = r_XMITloc;
    
    h = NaN(nh,nf);
    ts = NaN(nf,1);
    te = NaN(nf,1);
    par = NaN(nh,4,nf);
    parstd = NaN(nh,4,nf);
    azel = NaN(nf,2);
    
    for k=1:nf
        
        if ~mergedfile
            load(fp{k});
        else
            fname = ['data' flist(k).fname((end-11):(end-4))];
            tmpdata = getfield(load(fitdir,fname),fname);
            filefields = fieldnames(tmpdata);
            for ii=1:length(filefields)
                eval([char(filefields(ii)) '=tmpdata.' char(filefields(ii)) ';']);
            end
        end
        
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
            % if this created duplicates, just fill from the bottom
            if length(ih) ~= length(unique(ih))
                ih = 1:nh1;
            end
        end
        
        h(ih,k) = r_h;
        par(ih,:,k) = r_param(:,[1 2 3 5]);
        parstd(ih,:,k) = r_error(:,[1 2 3 5]);
        
        % remove failed iterations
        %    failed = (r_status ~= 0) | ( r_res(:,1) > 10 );
        %    failed = (r_status ~= 0) | ( r_res(:,1) > 30 );
        failed = ((r_status ~= 0) & (r_status ~=3)) | ( r_res(:,1) > 30 );
        par(failed,:,k) = NaN;
        parstd(failed,:,k) = NaN;
        
        azel(k,:) = [r_az r_el];
    end
    
    par(:,3,:) = par(:,3,:) .* par(:,2,:);
    parstd(:,3,:) = parstd(:,2,:).*par(:,3,:) + parstd(:,3,:).*par(:,2,:);
    
    
    % remove other than field-aligned data
    %
    % local field-aligned direction in E region
    [~,~, D, I,~,~,~,~,~,~] = igrfmagm(110000, loc(1), loc(2), decyear(r_time(1,:)));
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
    %    rminds = abs(mod(azel(:,1),360) - 187) > 3 | abs(abs(90-azel(:,2))-12.45) > 3;
    %    rminds = abs(mod(azel(:,1),360) - FAaz) > abs(FAdev) | abs(abs(90-azel(:,2))-(90-FAele)) > abs(FAdev);
    FAvec = [cos(FAele*pi/180)*cos(FAaz*pi/180) cos(FAele*pi/180)*sin(FAaz*pi/180) sin(FAele*pi/180)];
    BEAMvecs = [cos(azel(:,2)*pi/180).*cos(azel(:,1)*pi/180) cos(azel(:,2)*pi/180).*sin(azel(:,1)*pi/180) sin(azel(:,2)*pi/180)];
    rminds = acos(BEAMvecs*FAvec')*180/pi > FAdev;
    h(:,rminds) = [];
    ts(rminds) = [];
    te(rminds) = [];
    par(:,:,rminds) = [];
    parstd(:,:,rminds) = [];
    azel(rminds,:) = [];
    
    
    % do not accept Ti > 300 K below 100 km altitude or Ti < 50 K anywhere. If Ne < 1e9, keep the electron densities but do not use the temperatures. remove negative temperatures...
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
% list possible merged mat files and add their contents
fmerged = dir(fullfile(fitdir,'GUISDAP*.mat'));
if ~isempty(fmerged)
    for k=1:length(fmerged)
        disp(fullfile(fmerged(k).folder,fmerged(k).name))
        [hm,tsm,tem,parm,parstdm,locm,azelm,Im] = readGUISDAPpar(fullfile(fmerged(k).folder,fmerged(k).name),FAdev);
        h = [h,hm];
        ts = [ts;tsm];
        te = [te;tem];
        par = cat(3,par,parm);
        parstd = cat(3,parstd,parstdm);
        % loc and locm are identical if found. 
        if isempty(loc)
            loc = locm;
        end
        azel = [azel;azelm];
        % I and Im are identical if found.
        if isempty(I)
            I = Im;
        end
        
    end
    [~,iunique] = unique(ts);
    h = h(:,iunique);
    ts = ts(iunique);
    te = te(iunique);
    par = par(:,:,iunique);
    parstd = parstd(:,:,iunique);

    [~,isort] = sort(ts);
    h = h(:,isort);
    ts = ts(isort);
    te = te(isort);
    par = par(:,:,isort);
    parstd = parstd(:,:,isort);
end

end

