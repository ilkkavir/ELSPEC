function simudata = makeTestData(polycoefs,E,t,h,stdnoise)
%
% Create test data using fluxes corresponding the polynomial
% coefficients in polycoefs. parameters par are taken from the IRI
% model at times t and heights h. normally distributed pseudorandom
% noise with standard deviation stdnoise is added to the electron
% density profiles
%
% out = makeTestData(polycoefs,t)
%
%
%
%
%
%
%
%
%
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

nt = length(t);
nE = length(E)-1;
nh = length(h);

Ie = NaN(nE,nt);
neEnd = NaN(nh,nt);
ne = NaN(nh,nt);
dt = diff(t);
dt = [dt(1);dt(:)];
par = NaN(nh,4,nt);
parstd = NaN(nh,4,nt);
ts = t(:);
te = t(:)+dt;

% model parameters
% collect the model data in an array
model = NaN(nh,10,nt);
for it=1:nt
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = modelParams( t(it) , h );
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
end

% copy the model data to the par array
par(:,1,:) = model(:,8,:) + model(:,9,:) + model(:,10,:);
par(:,2,:) = model(:,2,:);
par(:,3,:) = model(:,3,:);
par(:,4,:) = 0;


ne0 = zeros(nh,1);

for it = 1:nt

    % ion production
    [A,Ec,dE] = ion_production( E(:) , h'*1000 , model(:,4,it) , ...
                                model(:,5,it) , model(:,6,it) , model(:,7,it) , ...
                                model(:,1,it) , 'Fang' );
    A(isnan(A)) = 0;

    % effective recombination rate
    alpha = effective_recombination_coefficient( h , ...
                         model(:,3,it) , ...
                         model(:,9,it) , ...
                         model(:,9,it).*0 , ...
                         model(:,8,it) , ...
                         'Rees' );

    % calculate the model spectra
    Ie(:,it) = model_spectrum( polycoefs( : , it ) , Ec );

    % ne at end of integration
    neEnd(:,it) = integrate_continuity_equation( ne0 , Ie(:,it) , ...
                                                 dt(it) , A , dE(:) , ...
                                                 alpha , 'endNe');
    % average ne
    ne(:,it) = integrate_continuity_equation( ne0 , Ie(:,it) ...
                                                      , dt(it) , A ...
                                                      , dE(:) ...
                                                      , alpha ...
                                                      , 'integrate' ...
                                                      );
    %                                                      , 'equilibrium' ...
    ne0 = ne(:,it);

end



rng('shuffle');
pp = ne + stdnoise*randn(size(ne));
ppstd = zeros(size(pp)) + stdnoise;

pp(isnan(pp)) = 0;
ppstd(isnan(ppstd)) = 1e15;

simudata.h = h';
simudata.ts = ts;
simudata.te = te;
simudata.pp = pp;
simudata.ppstd = ppstd;
simudata.par = par;
simudata.parstd = parstd;
simudata.iri = model;
simudata.Ie = Ie;

end


