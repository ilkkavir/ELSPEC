function [h,ts,te,pp,ppstd,par,parstd,model,f107,f107a,f107p,ap,loc,azel,I] = readFitData( ppdir , ...
                                                  fitdir , hmin , ...
                                                  hmax , tmin , tmax ...
                                                  , exp  , radar , ...
                                                  version , tres , ...
                                                      readIRI , FAdev , BottomStdFail)
%
%  Read GUISDAP raw densities (power profiles), GUISDAP fit
%  results, and model (IRI and MSIS) parameters.
%
% [h,ts,te,pp,ppstd,par,parstd,model,f107,f107a,f107p,ap,loc,azel.I] =
%    readFitData( ppdir , fitdir , hmin , hmax , tmin , tmax  , exp
%    , radar , version , tres , readIRI , FAdev)
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
%  radar   EISCAT radar name ['uhf','vhf','esr']
%  version EISCAT experiment version number [1,2,3,...]
%  tres    "type" of time resolution 'best' or 'dump'
%  readIRI logical, do we read the tabulated IRI data at all
%  FAdev   maximum beam direction deviation from field-aligned  [deg], default 3
%  BottomStdFail Standard deviation given for the model Ne in the lowest observed gate when
%                the GUISDAP fit has failed. 
%
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
%  model   length(h) x 10 x length(ts) array of model parameters
%  f107
%  f107a
%  f107p
%  ap
%  loc     radar coordinates (lat, lon) [deg]
%  azel    radar beam azimuths and elevations [deg]
%  I       magnetic inclination [deg]
%
% The four plasma parameters are Ne [m^-3],Ti [K], Te [k], and Vi
% [ms^-1]. Failed iterations and results with chi-squared larger
% than 10 are replaced with NaNs.
%
% The 10 model parameters are Tn [K],Ti [K] ,Te [K], nN2 [m^-3],
% nO2 [m^-3], nO [m^-3], nAr[m^-3], nNOp [m^-3], nO2p [m^-3], nOp [m^-3],
% where 'n' refers to number density and 'p' to a positive ion.
%
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

if (hmin<80)
    error('Minimum height must be at least 80 km.');
end
if (hmax>150)
    error('Maximum height must be at most 150 km.');
end

% read GUISDAP power profiles and fit results
[h,ts,te,pp,ppstd,par,parstd,loc,azel,I] = readGUISDAPdata( ppdir , ...
                                                  fitdir , hmin , ...
                                                  hmax , tmin , tmax ...
                                                  , exp , radar , ...
                                                 version , tres , FAdev );

% collect the model data in an array
model = NaN(length(h),10,length(ts));
f107 = NaN(length(ts),1);
f107a = NaN(length(ts),1);
f107p = NaN(length(ts),1);
ap = NaN(length(ts),1);

% no need to continue if there was no data
if isempty(ts) | isempty(h)
    return
end

for it=1:length(ts)
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp,f107it,f107ait,apit] = modelParams2( ...
        (ts(it)+te(it))./2 , h , loc , readIRI );
    [~,~,~,~,~,~,~,~,~,~,f107pit,~,~] = modelParams( ...
        (ts(it)+te(it))./2-86400 , h , loc , readIRI );
    model(:,1,it) = Tn;
    model(:,2,it) = Ti;
    model(:,3,it) = Te;
    model(:,4,it) = nN2;
    model(:,5,it) = nO2;
    model(:,6,it) = nO;
    model(:,7,it) = nAr;
    model(:,8,it) = nNOp;
    model(:,9,it) = nO2p;
    model(:,10,it) = nOp;

    f107(it) = f107it;
    f107a(it) = f107ait;
    f107p(it) = f107pit;
    ap(it) = apit;
end

% replace NaN values with interpolation or model and large errors

% find the NaN points, notice the order of the incides...
[i1,i3,i2] = ind2sub(size(par),find(isnan(par)));

% first try to fill with linear interpolation in height, but excluding end points
par = fillmissing(par,'linear',1,'endvalues',NaN);

% then fill the remaining NaNs with model values. Give a large
% variance to both the model values and the interpolated ones
% Ne
%[i1,i2] = find(isnan(par(:,1,:)));

% replace the remaining NaN's with model values
ii1 = i1(i3==1);
ii2 = i2(i3==1);
for inan = 1:length(ii1)
    if isnan(par(ii1(inan),1,ii2(inan)))
        if readIRI
            par(ii1(inan),1,ii2(inan)) = model(ii1(inan),8,ii2(inan)) + ...
                model(ii1(inan),9,ii2(inan)) + model(ii1(inan),10, ...
                                                     ii2(inan));
        else
            par(ii1(inan),1,ii2(inan)) = 0;
        end
    end
    parstd(ii1(inan),1,ii2(inan)) = 1e12;
end

% Ti
%[i1,i2] = find(isnan(par(:,2,:)));
ii1 = i1(i3==2);
ii2 = i2(i3==2);
for inan = 1:length(ii1)
    if isnan(par(ii1(inan),2,ii2(inan)))
        par(ii1(inan),2,ii2(inan)) = model(ii1(inan),2,ii2(inan));
    end
    parstd(ii1(inan),2,ii2(inan)) = 1e4;
end
% Te
%[i1,i2] = find(isnan(par(:,3,:)));
ii1 = i1(i3==3);
ii2 = i2(i3==3);
for inan = 1:length(ii1)
    if isnan(par(ii1(inan),3,ii2(inan)))
        par(ii1(inan),3,ii2(inan)) = model(ii1(inan),3,ii2(inan));
    end
    parstd(ii1(inan),3,ii2(inan)) = 1e4;
end
% Vi
%[i1,i2] = find(isnan(par(:,4,:)));
ii1 = i1(i3==4);
ii2 = i2(i3==4);
for inan = 1:length(ii1)
    if isnan(par(ii1(inan),4,ii2(inan)))
        par(ii1(inan),4,ii2(inan)) = 0;
    end
    parstd(ii1(inan),4,ii2(inan)) = 1e4;
end

% pp
[i1,i2] = find(isnan(pp));
pp = fillmissing(pp,'linear',1,'endvalues',NaN);
for inan = 1:length(i1)
    if isnan(pp(i1(inan),i2(inan)))
        if readIRI
            pp(i1(inan),i2(inan)) = model(i1(inan),8,i2(inan)) + ...
                model(i1(inan),9,i2(inan)) + model(i1(inan),10,i2(inan));
        else        
            pp(i1(inan),i2(inan)) = 0;
        end
    end
    ppstd(i1(inan),i2(inan)) = 1e12;
end

% remove zero variances in pp...
izerostd = ppstd<=0;
ppstd(izerostd) = 1e12;


% If variance 1e12 was put to the lowest gate, assume that Ne must have been small and put
% there a smaller variance to avoid unrealistic Ne profiles when the bottom part is very noisy
% do this in the lowest gate from which we have at least one measurement
ilow = 1;
while all(ppstd(ilow,:)==1e12)
    ilow = ilow + 1;
end

ppstd(ilow,ppstd(ilow,:)==1e12) = BottomStdFail;



end
