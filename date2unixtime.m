function [unix_time] = date2unixtime(utc_date)
%
% date2unixtime convert utc date into unix time 
% (seconds since 1970-01-01 00:00:00 UTC)
%
% unix_time = date2unixtime(utc_date)
%
% IV 2016
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

unix0 = datenum(1970,01,01); %Start time of unix time as matlab datenum
unixday = datenum(utc_date(1),utc_date(2),utc_date(3)); % matlab datenum of utc_date

unix_time = 60 * ( 60 * ( 24 * ( unixday - unix0 ) + utc_date(4) ) ...
                   + utc_date(5) ) + utc_date(6);

end
