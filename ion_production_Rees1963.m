function [q,Ec,dE] = ion_production_Rees1963(E,z,nN2,nO2,nO,nAr,Tn)
%   Ionization profiles of monoenergetic electron precipitation
%   according to Rees [1963].
%
%   [q,Ec,dE] = ion_production_Rees1963(E,z,nN2,nO2,nO,Tn)
%
% INPUT:
%   E    electron energies [eV] (1 x (nE+1) array)
%   z    altitudes [m] (nz x 1 array)
%   nN2  molecular Nitrogen density [m^-3] (nz x 1 array)
%   nO2  molecular Oxygen density [m^-3] (nz x 1 array)
%   nO   atomic oxygen density [m^-3] (nz x 1 array)
%   nAr   Argon density [m^-3] (nz x 1 array)
%   Tn   neutral temperature [K] (nz x 1 array)
%
% OUTPUT:
%   q   ionization profiles [m^-1] (nz x nE array)
%       at centres of the energy bins, i.e. energies (E(1)+E(2))/2,
%       (E(2)+E(3))/2 etc.
%   Ec  centre points of the energy bins [eV]
%   dE  widths of the energy bins [eV]
%
% I. Virtanen 2018
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

% number of height gates
nz = length(z);

% number of energy bins. One less than input points, because our steps are from E(i) to E(i+1), etc.
nE = length(E) - 1;

% widths of energy bins
dE = E(2:end) - E(1:nE);

% centre points of the energy bins
Ec = E(1:nE) + dE/2;

% model atmosphere from Rees' table 1. Numbers copied from CARD
% norlib_cray.f, function PRODU
AA = [60.,64.,68.,72.,76.,80.,84.,88.,92.,96.,100.,108.,114.,120., ...
      126.,132.,140.,150.,162.,192.,230.,252.,276.,300.];

BBZ = [2.34E-1,1.37E-1,7.76E-2,4.20E-2,2.17E-2,1.08E-2,5.16E-3, ...
       2.46E-3,1.19E-3,5.88E-4,3.10E-4,1.04E-4,5.18E-5,2.75E-5, ...
       1.65E-5,1.08E-5,6.90E-6,4.40E-6,2.96E-6,1.41E-6,7.05E-7, ...
       5.02E-7,3.59E-7,2.65E-7];


BBN = [6.33E15,3.93E15,3.39E15,1.39E15,7.72E14,4.04E14,1.99E14, ...
       9.48E13,4.37E13,2.07E13,1.04E13,3.18E12,1.43E12,6.61E11, ...
       3.40E11,1.91E11,9.70E10,4.92E10,2.66E10,9.50E09,3.89E09, ...
       2.57E09,1.72E09,1.20E09];

XXL1 = [.00,.04,.08,.14,.24,.42,.65,1.00,1.63,1.91,1.85,1.70,1.52, ...
        1.33,1.17,1.05,.95,.84,.73,.62,.52,.42,.33,.25,.17,.12,.06, ...
        .02,.00,.00];

RATI = -.4 + .05*(0:29); % x-axis of Table 1 in Rees 1963.

% Average energy per ion-electron-pair [keV]
dEion = .035;


m_e     = 9.10938291e-31;           % electron rest mass [kg]
m_p     = 1.672621778e-27;          % proton rest mass [kg]
m_n     = 1.674927352e-27;          % neutron rest mass [kg]
% mass density (g/cm^3)
ro = ((nO*8*(m_n+m_p+m_e) + nO2*2*8*(m_n+m_p+m_e) + nN2*2*7.5*(m_n+m_p+m_e)) ...
     + nAr*(18*m_p+22*m_n+18*m_e))*1e-3;

q = zeros(nz,nE);

for iE = 1:nE
    E0 = Ec(iE)/1e3;                         % Rescale electron energy to keV
    R = 4.57e-6 * E0^1.75;                   % penetration depth
                                             % g/cm^2
    %R0 = R./ro';                             % this is copied from
                                             % CARD, but I think it
                                             % may be wrong...
    zz = interp1(AA,BBZ,z/1000,'spline',0);  % atmospheric depths (g/cm^2)

    znz = interp1(BBZ,BBN,zz,'spline',0);    % number densities at z
                                             % (cm^-3)
    znr = interp1(BBZ,BBN,R,'spline',0);     % number density at
                                             % R (cm^-3)
    zrr = interp1(zz,ro,R,'spline');                  % mass density at R

    ratn = znz./znr;                         % number density ratio
                                             % on eq 1 of Rees
                                             % 1963.
    ratz = zz./R;                            % ratio of atmospheric depths

    R0 = R./zrr;

    % normalized energy dissipation distribution function
    lambda = interp1(RATI,XXL1,ratz,'spline',0);

    % the production rates.
    q(:,iE) = (E0./R0).*(lambda./dEion).*ratn.*100;

end
%    q = q.*2;

q(q<0) = 0;
q(isnan(q)) = 0;

end