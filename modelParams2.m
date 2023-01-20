function [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp,f107,f107a,ap] = modelParams2( time , heights , loc , readIRI )
%
% IRI and MSIS model parameters. The parameters are read from a
% table delivered as part of the ElSpec package.
%
% [Tn,Ti,Te,nN2,nO2,nO,nAr,nO2p,nN2p,nNOp] = modelParams( time ,
% heights , loc , readIRI )
%
% INPUT
%   time     time as unix time
%   heights  heights in km
%   loc      [lat lon] of the radar site
%   readIRI  logical, do we read the tabulated IRI data?
%
%
% OUTPUT:
%  Tn   neutral temperature
%  Ti   ion temperature
%  Te   electron temperature
%  nN2  N2 number density
%  nO2  O2 number density
%  nO   O number density
%  nAr  Ar number density
%  nNOp NO+ ion number density
%  nO2p O2+ ion number density
%  nOp  O+ ion number density
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

persistent year_prev iripar apf107table readIRIprev locprev

%% time rounded to the closest full hour
%%matlabtime = datetime(round(time/3600)*3600,'ConvertFrom','posixtime');

% time as matlab datetime. The tabulated IRI data are at half hours, so no need to round anymore
matlabtime = datetime(round(time/3600)*3600,'ConvertFrom','posixtime');
year = matlabtime.Year;
month = matlabtime.Month;
hour = matlabtime.Hour;
%if hour==0
%    hour=1; % FIX THIS!!
%end


if readIRI

    if isempty(year_prev) | year_prev ~= year | readIRIprev ~= readIRI | locprev ~= loc
        if sum(abs(loc(1:2) - [69.58 19.23])) < 5
            load(['ElSpecModelParametersTRO',num2str(year),'.mat']);
        elseif sum(abs(loc(1:2) - [78.15 16.02])) < 5
            load(['ElSpecModelParametersESR',num2str(year),'.mat']);
        else
            error(['IRI parameters are not available for location ' loc])
        end
    end


    %hh = round((heights-80)/2)+1;

    if (any(heights<80 | heights>150))
        error('Heights must be between 80 and 150 km');
    end

    partmp = squeeze(iripar(month,matlabtime.Day,hour+1,:,:));

    hh = 80:2:150;

    nOp = interp1(hh,partmp(:,1),heights);
    nO2p = interp1(hh,partmp(:,2),heights);
    nNOp = interp1(hh,partmp(:,3),heights);
    Ti = interp1(hh,partmp(:,4),heights);
    Te = interp1(hh,partmp(:,5),heights);

else
    nOp = heights*NaN;
    nO2p = heights*NaN;
    nNOp = heights*NaN;
    % Tn is copied to Ti and Te at end of this function.
end

% 20170508
% use matlab atmosnrlmsise00 instead of the nrlmsise00 within IRI
% for neutrals

% first read the ap-indices and f017 values from file apf107.dat
% the table apf107table is a persistent variable, so it is enough
% to read the original file once.
%
% The latest apf107.dat file is available from https://chain-new.chain-project.net/echaim_downloads/apf107.dat
%
%
% file format:
%
%  year(I3), month(I3), day(I3), 3-hour Ap indices for the UT
%  intervals (0-3), )3-6), )6-9), .., )18-21), )21-24(  in an array
%  of dimension 8 (8I3), daily Ap (I3), -11(I3),
%  F10.7 radio flux for the day (F5.1), 81-day average of F10.7
%  radio flux (F5.1), 365-day average of F10.7 centered on the date
%  of interest (F5.1). At start and end of the index file the
%  81-day and 365-day averages are calculated taking only the
%  available indices, e.g. for the first date the 81-day average is
%  only over 40 F10.7 values and over 41 values on the 2nd date.

if isempty(apf107table)
    apfline = 1;
    apf107fid = fopen('apf107.dat');
    while ~feof(apf107fid)
        apf107str = fgetl(apf107fid);
        %        apf107table(apfline,:)=sscanf(fgetl(apf107fid),
        %        '%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%5f%5f%5f');
        apf107table(apfline,1) = str2double(apf107str(1:3));
        apf107table(apfline,2) = str2double(apf107str(4:6));
        apf107table(apfline,3) = str2double(apf107str(7:9));
        apf107table(apfline,4) = str2double(apf107str(10:12));
        apf107table(apfline,5) = str2double(apf107str(13:15));
        apf107table(apfline,6) = str2double(apf107str(16:18));
        apf107table(apfline,7) = str2double(apf107str(19:21));
        apf107table(apfline,8) = str2double(apf107str(22:24));
        apf107table(apfline,9) = str2double(apf107str(25:27));
        apf107table(apfline,10) = str2double(apf107str(28:30));
        apf107table(apfline,11) = str2double(apf107str(31:33));
        apf107table(apfline,12) = str2double(apf107str(34:36));
        apf107table(apfline,13) = str2double(apf107str(37:39));
        apf107table(apfline,14) = str2double(apf107str(40:44));
        apf107table(apfline,15) = str2double(apf107str(45:49));
        apf107table(apfline,16) = str2double(apf107str(50:54));
        apfline = apfline+1;
    end
    fclose(apf107fid);
end

% the two-digit year number use in the apf107 file...
ynum = mod(matlabtime.Year,100);
% day-of-year
dayofyear = day(matlabtime,'dayofyear');
% second-of-day
%secofday = ((matlabtime.Hour*60) + matlabtime.Minute*60) + matlabtime.Second;
secofday = seconds(matlabtime-dateshift(matlabtime,'start','day'));

% the correct line from apf107table
apf107line = find( apf107table(:,1)==ynum & apf107table(:,2)== ...
                   matlabtime.Month & apf107table(:,3)==matlabtime.Day);

if ~isempty(apf107line)
% pick the f10.7 values
f107 = apf107table(apf107line,14);
f107a = apf107table(apf107line,15);

% daily ap
ap = apf107table(apf107line,12);

% pick the relevant part of the 3-hour ap indices
ap3 = apf107table(apf107line-3:apf107line,4:11)'; % everything we
                                                  % might need
ap3 = ap3(:); % convert into a vector (note the "'" in previous line)
ap3 = ap3(1:(end-(8-ceil(matlabtime.Hour/3)))); % remove extra from
                                                % the end
ap3 = ap3( end-19 : end); % remove extra from the beginning

% then form the aph input for atmosnrlmsise00
aph = [ap,ap3(end),ap3(end-1),ap3(end-2),ap3(end-3),mean(ap3(end-4: ...
                                                  end-11)),mean(ap3(end-12:end-19))];

%[Tmsis,rhomsis] = atmosnrlmsise00( heights*1000,69.58646 , 19.22743, ...
%                                   matlabtime.Year, dayofyear, secofday, f107a, f107,aph);
[Tmsis,rhomsis] = atmosnrlmsise00( heights*1000,loc(1) , loc(2), ...
                                   matlabtime.Year, dayofyear, secofday, f107a, f107,aph);

else
    %[Tmsis,rhomsis] = atmosnrlmsise00( heights*1000,69.58646 , 19.22743, ...
    %                                   matlabtime.Year, dayofyear, secofday);
[Tmsis,rhomsis] = atmosnrlmsise00( heights*1000,loc(1) , loc(2), ...
                                   matlabtime.Year, dayofyear, secofday);
end

nO = rhomsis(:,2);
nN2 = rhomsis(:,3);
nO2 = rhomsis(:,4);
nAr = rhomsis(:,5);
Tn = Tmsis(:,2);

% copy Tn to Ti and Te if IRI was not used
if ~readIRI
    Te = Tn;
    Ti = Tn;
end

year_prev = year;
readIRIprev = readIRI;
locprev = loc;

end