function [q,qnorm,Ec,dE,Emax] = calculateIonProduction(time,heights,loc,Ec,model)
%
% [q,Ec,dE,Emax] = calculateIonProduction(time,heights,loc,Ec,model)
%
% Calculate the ion production by mono-energetic electrons based on time and location.
%
% INPUTS:
%  time    time as matlab datetime
%  heights heights in km
%  loc     latitude and longitude (geodetic) in degrees
%  Ec      centre points of the energy bins in eV
%  model   "fang", "sergienko", or "rees"
%
%
% OUPUTS:
%  q      the ion production rates
%  Ec     the input energy grid
%  dE     widths of the energy bins
%  Emax   energy with highest production at each height
%
%
%
% EISCAT site coordinates (geographic):
% TRO: [69.58 19.23]
% KIR: [67.87 20.43]
% SOS: [67.37 26.62]
% ESR: [78.15 16.02]
%
%
% IV 2022
%
%

% neutral atmosphere from msis
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp,f107,f107a,ap] = modelParams( posixtime(time) , heights , loc , 0);

    % magnetic inclination
    [~,~, ~, I,~,~,~,~,~,~] = igrfmagm(110000, loc(1), loc(2), decyear(time));

    % the production rates
    [q,Ec,dE] = ion_production(Ec,heights*1000,nN2,nO2,nO,nAr,Tn,model,I);

    qnorm = q;
    for k=1:length(Ec)
        qnorm(:,k) = q(:,k)/max(q(:,k));
    end
    
    % find the energy that has highest prodution at this altitude
    Emax = [];
    for k=1:length(heights)
        [~,imax] = max(qnorm(k,:));
        Emax(k) = Ec(imax);
    end

end
    