function [prod , hh] = modelProductionProfile(Ec,time,lat,lon,ionomodel)
%
% Calculate ion production rate as function of altitude at 80-150 km altitudes
% at the given time and location
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later
%


    % altitude grid
    hh= 80:.1:150;

    % time as unix time
    tt = posixtime(datetime(time));
    
    % call the models
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp,f107,f107a,ap] = modelParams2( ...
        tt , hh , [lat,lon,0] , false );
    [~,~,~,~,~,~,~,~,~,~,f107p,~,~] = modelParams2( ...
        tt-86400 , hh , [lat,lon,0] , false );

    [~,~, D, I,~,~,~,~,~,~] = igrfmagm(110000, lat, lon, decyear(time));

    prod = ion_production([-1 1]+Ec,hh,nN2,nO2,nO,nAr,Tn,ionomodel,I);
    
    
end
