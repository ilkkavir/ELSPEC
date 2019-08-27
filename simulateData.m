function simudata = simulateData(nestd,ne0std,Sstd)
%
% Create simulated electron density data.
%
% Assumes a Maxwellian spectrum, whose peak energy and flux as
% function of time are hard-coded in the function
%
%
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


% one hour, 5s steps
t = (0:720)*5;

% start and end times, assume start at 20151109T2000
ts = t + 1447099200;
te = ts + 5;

% a 2 keV peak centred ad t=800s
ep2 = -(((t-800)/10).^2)+2;
ep2(ep2<0)=0;

% 5 keV centred at t=1400s
ep5 = -(((t-1400)/6).^2)+5;
ep5(ep5<0)=0;

% 10 keV centred at t=2000s
ep10 = -(((t-2000)/7).^2)+10;
ep10(ep10<0)=0;


% 15 keV centred at t=2600s
ep15 = -(((t-2600)/12).^2)+15;
ep15(ep15<0)=0;

% 20 keV centred at t=3200s
ep20 = -(((t-3200)/3).^2)+20;
ep20(ep20<0)=0;


% the combined characteristic energy, .5 keV background
ep = ep2 + ep5 + ep10 + ep15 + ep20;
ep = ep + .5;
ep = ep * 1000; % conversion to keV

% an oscillating total flux (eVm^-2s^-1)
Q = (sin(t/100)+1.01)*.005 / 1.60217662e-19;


% energy grid in eV
E = logspace(log10(10),log10(1e6),500);

% heights
h = 80:1.5:150;

nh = length(h);
nt = length(ts);

% the model parameters
% collect the model data in an array
model = NaN(length(h),10,length(ts));
for it=1:length(ts)
    [Tn,Ti,Te,nN2,nO2,nO,nAr,nNOp,nO2p,nOp] = modelParams( ...
        (ts(it)+te(it))./2 , h );
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

pp = NaN(nh,nt);
nE = length(E)-1;
simudata.Ie = NaN(nE,nt);
simudata.alpha = NaN(nh,nt);
simudata.q = NaN(nh,nt);
for tt=1:nt
    % ion production
    [A,Ec,dE] = ion_production( E , h*1000 , model(:,4,tt) , ...
                                model(:,5,tt) , model(:,6,tt) , model(:,7,tt) , ...
                                model(:,1,tt) , 'Fang' );
    A(isnan(A)) = 0;


    % recombination
    alpha = effective_recombination_coefficient( h , ...
                                                 model(:,3,tt) , ...
                                                 model(:,9,tt) , ...
                                                 model(:,9,tt).*0 , ...
                                                 model(:,8,tt) , ...
                                                 'Rees' );

    % the maxwellian flux
    S = MaxwellFlux(ep(tt),Q(tt),Ec);%Q(tt)/(2*ep(tt).^3)*Ec.*exp(-Ec/ep(tt));

    % a kappa-distribution...
    %    S = KappaFlux(ep(tt),Q(tt),5,Ec);%Q(tt)/(2*ep(tt).^3)*(kappa-1)*(kappa-2)/kappa^2*Ec*(1+Ec/(kappa*ep(tt)));

    % add some noise to the spectrum
    S = S + randn(size(S)).*S*Sstd;
    S(S<0) = 0;

    if tt==1
        ne0 = zeros(nh,1);
    end
    pp(:,tt) = integrate_continuity_equation(ne0(:),S(:),5,A,dE(:), ...
                                             alpha(:),'integrate');
    ne0 = integrate_continuity_equation(ne0(:),S(:),5,A,dE(:),alpha(:), ...
                                        'endNe') + randn(size(ne0))*ne0std;
    ne0(ne0<1) = 1;

    simudata.Ie(:,tt) = S;
    simudata.alpha(:,tt) = alpha;
    simudata.q(:,tt) = A*(simudata.Ie(:,tt).*dE');

end

% add noise to pp
pp = pp + randn(size(pp))*nestd;

simudata.h = h;
simudata.ts = ts;
simudata.te = te;
simudata.pp = pp;
simudata.ppstd = ones(nh,nt)*nestd;

simudata.par = model(:,1:4,:);
simudata.par(:,1,:) = pp; % ne
simudata.par(:,4,:) = pp.*0; % vi (not used...)
simudata.parstd = simudata.par;
simudata.parstd(:,1,:) = nestd;
simudata.parstd(:,2:4,:) = NaN;

simudata.iri = model;




simudata.hmin = 80;
simudata.hmax = 150;
simudata.iono_model = 'Fang';
simudata.recomb_model = 'Rees';
simudata.integ_type = 'integrate';
simudata.E = E;
simudata.Emin = 1e3;
simudata.ninteg = 1;
simudata.nstep = 1;
simudata.Ec = Ec;
simudata.dE = dE;
simudata.dt = diff(simudata.te);
simudata.dt = [simudata.dt(1);simudata.dt(:)];
simudata.ne = simudata.pp;
simudata.nestd = simudata.ppstd;
Eind_fac = simudata.Ec >= simudata.Emin;
EdE = simudata.dE(Eind_fac)'.*simudata.Ec(Eind_fac)';
for tt=1:nt
    simudata.FAC(tt) = sum(simudata.Ie(Eind_fac,tt).* ...
                           simudata.dE(Eind_fac)')*1.60217662e-19;
    simudata.Pe(tt) = sum(simudata.Ie(Eind_fac,tt).*EdE)* ...
        1.60217662e-19;

end
simudata.FACstd = simudata.FAC*0;
simudata.PeStd = simudata.Pe*0;
simudata.A = A;

simudata.chisqr = simudata.Pe*NaN;
