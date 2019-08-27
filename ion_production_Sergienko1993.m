function [q,Ec,dE] = ion_production_Sergienko1993(E,z,nN2,nO2,nO,Tn)
%
%   Ionization profiles of monoenergetic electron precipitation
%   according to Sergienko and Ivanov [1993]. This function calculates the
%   ionization profiles for electron precipitation into the Earths
%   atmosphere.
%
% CALLING:
%  [q,Ec,dE] = ion_production_Sergienko1993(E,z,nN2,nO2,nO,Tn)
%  q = ion_production_Sergienko1993(E,z,nN2,nO2,nO,Tn)
%
% INPUT:
%   E    electron energies [eV] (1 x nE array)
%   z    altitudes [m] (nz x 1 array)
%   nN2  molecular Nitrogen density [m^-3] (nz x 1 array)
%   nO2  molecular Oxygen density [m^-3] (nz x 1 array)
%   nO   atomic oxygen density [m^-3] (nz x 1 array)
%   Tn   neutral temperature [K] (nz x 1 array)
%
% OUTPUT:
%   q   ionization profiles [m^-1] (nz x nE array)
%       at centres of the energy bins, i.e. energies (E(1)+E(2))/2,
%       (E(2)+E(3))/2 etc.
%   Ec  centre points of the energy bins [eV]
%   dE  widths of the energy bins [eV]
%
%   Based on the function ionization_profiles_from_flux by B. Gustavsson.
%
% I. Virtanen 2016
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% number of energy bins. One less than input points, because our steps are from E(i) to E(i+1), etc.
nE = length(E) - 1;

% widths of energy bins
dE = E(2:end) - E(1:nE);

% centre points of the energy bins
Ec = E(1:nE) + dE/2;

% number of height points
nz = length(z);

% Physical constants
kB	= 1.380662e-23;		    % Boltzmann constant [J/K]
Re	= 6.378e6;		    % Radius of earth [m]
m_e     = 9.10938291e-31;           % electron rest mass [kg]
m_p     = 1.672621778e-27;          % proton rest mass [kg]
m_n     = 1.674927352e-27;          % neutron rest mass [kg]

% mass density (kg/m^3)
rho = (nO*8*(m_n+m_p+m_e) + nO2*2*8*(m_n+m_p+m_e) + nN2*2*7.5*(m_n+m_p+m_e));

q = zeros(nz,nE);

% energy loss per ion pair for N2, O2, and O
% Sergienko & Ivanov, table 1
E_cost_ion = [36.8,28.2,26.8];

% from Monte-Carlo method of Sergienko & Ivanov (eq. 7)
% kN2, kO2, kO
ki = [1 0.7 0.4];

% Pi in eq 7 of Sergienko & Ivanov
Partitioning = ( [ki(1)*nN2,ki(2)*nO2,ki(3)*nO] ./ ...
                 repmat(sum([ki(1)*nN2,ki(2)*nO2,ki(3)*nO],2),1,3) );

% First calculate the energy deposition as a function of altitude
for iE = 1:nE,

  % for isotropic electron velocity distributions
  A = energy_degradation_Sergienko1993(Ec(iE),z,nN2,nO2,nO,rho)';

  % ion production rate by precipitating electrons [m^-1]
  % q will need to be multplied with electron flux [(m^2 eV s)^-1]
  % and energy bin widths [eV] to get the ion production rate
  % [(m^3 s)^-1]
  q(:,iE) = A.*(Partitioning(:,1)/E_cost_ion(1) + ...
                Partitioning(:,2)/E_cost_ion(2) + ...
                Partitioning(:,3)/E_cost_ion(3));

end

end % ion_production_Sergienko1993


function [A] = energy_degradation_Sergienko1993(E,z,nN2,nO2,nO,rho);
%
% Energy degradation of electrons according to Sergienko et al. (1993).
%
% based on the function energy_deg (by B. Gustavsson?)
%
%  A = energy_degradation_Sergienko1993(Energy,z,nN2,nO2,nO,density);
%
%
% INPUT:
%   E        electron energy [eV]
%   z        altitudes [m] (nz x 1 array)
%   nN2      molecular Nitrogen density [m^-3] (nz x 1 array)
%   nO2      molecular Oxygen density [m^-3] (nz x 1 array)
%   nO       atomic oxygen density [m^-3] (nz x 1 array)
%   rho      mass density [kg/m^3] (nz x 1 array)
%
% OUTPUT:
%   A       energy degradation as function of altitude [eV m^-1]
%


nz = length(z);
zetm = zeros(nz,1);
dH = gradient(z);
for ih = (nz-1):-1:1,
  dzetm = ( rho(ih+1) + rho(ih) ) * dH(ih) / 2;
  zetm(ih) = zetm(ih+1)+dzetm; % [kg m^-2]
end

A = zeros(1,nz);

% Integral average range of the electron flux [kg m^-2]
% for isotropic velocity distributions. Eq. A1 and A3 of
% Sergienko & Ivanov
pr = 10*(1.64e-06*(E/1000).^1.67.*(1.+9.48e-02*(E/1000).^(-1.57)));
% a monodirectional beam would be
%pr = 10*(2.16e-06*(E/1000).^1.67.*(1.+9.48e-02*(E/1000).^(-1.57)));

for ih = 1:nz,

    % normalized distance from the electron source (dimensionless)
    hi = zetm(ih)/pr;

    % See equation A4 of Sergienko & Ivanov. The remaining terms
    % are added in ion_production_Sergienko1993
    A(ih) = rho(ih)*energy_dissipation_Sergienko1993(hi,E,2)*E*(1-albedo_Sergienko1993(E,2))/pr;

end

end % energy_degradation_Sergienko1993




function [alb]=albedo_Sergienko1993(Energy,pad),
%
% Albedo-flux of Sergienko & Ivanov
%
% these numbers are not tabulated in the paper, but they seem to
% match with those in Fig13b
%
logE_p=[1.69 1.8:0.1:3.7];

Param(1,:)=[0.352 0.344 0.334 0.320 0.300 0.280 0.260 0.238 0.218 0.198 0.180 0.160 0.143 0.127 0.119 0.113 0.108 0.104 0.102 0.101 0.100];
Param(2,:)=[0.500 0.492 0.484 0.473 0.463 0.453 0.443 0.433 0.423 0.413 0.403 0.395 0.388 0.379 0.378 0.377 0.377 0.377 0.377 0.377 0.377];

logE=log10(Energy);

if logE > 3.7,
  alb=Param(pad,end);
else
  alb=interp1(logE_p,Param(pad,:),logE,'linear','extrap');
end

end % albedo_Sergienko1993



function [lam]=energy_dissipation_Sergienko1993(hi,Energy,pad)

% see table 6 and equation A2 of Sergienko and Ivanov, A new approach to calculate
% the excitation of atmospheric gases by auroral electron impact
% annales geophysicae, 11, 717-727, 1993.

% monodirectional
logE_m=[1.69 1.8:0.1:3.7];
Param_m(1,:)=[1.43 1.51 1.58 1.62 1.51 1.54 1.18 1.02 0.85 0.69 0.52 0.35 0.21 0.104 0.065 0.05 0.04 0.03 0.03 0.025 0.021];
Param_m(2,:)=[0.83 0.77 0.72 0.67 0.63 0.59 0.56 0.525 0.495 0.465 0.44 0.42 0.40 0.386 0.37 0.36 0.35 0.34 0.335 0.325 0.32];
Param_m(3,:)=-[0.025 0.030 0.040 0.067 0.105 0.155 0.210 0.275 0.36 0.445 0.51 0.61 0.69 0.77 0.83 0.865 0.90 0.92 0.935 0.958 0.96];
Param_m(4,:)=[-1.67 -1.65 -1.62 -1.56 -1.46 -1.35 -1.20 -0.98 -0.70 -0.37 -0.063 0.39 0.62 0.92 1.11 1.25 1.36 1.44 1.50 1.55 1.56];

% isotropic
logE_i=[1.69 1.8:0.1:3.0];
Param_i(1,:)=[0.041 0.051 0.0615 0.071 0.081 0.09 0.099 0.1075 0.116 0.113 0.13 0.136 0.139 0.142];
Param_i(2,:)=[1.07 1.01 0.965 0.9 0.845 0.805 0.77 0.735 0.71 0.69 0.67 0.665 0.66 0.657];
Param_i(3,:)=-[0.064 0.1 0.132 0.171 0.2 0.221 0.238 0.252 0.261 0.267 0.271 0.274 0.276 0.277];
Param_i(4,:)=-[1.054 0.95 0.845 0.72 0.63 0.54 0.475 0.425 0.38 0.345 0.319 0.295 0.28 0.268];


if pad == 1,

    logE=log10(Energy);
    if Energy >= 5000,
        lam=(Param_m(1,end)*hi+Param_m(2,end)).*exp(Param_m(3,end).*hi.*hi+Param_m(4,end)*hi);
    else
        % equation A2 from Sergienko and Ivanov
        C1=interp1(logE_m,Param_m(1,:),logE);
        C2=interp1(logE_m,Param_m(2,:),logE);
        C3=interp1(logE_m,Param_m(3,:),logE);
        C4=interp1(logE_m,Param_m(4,:),logE);
        lam=(C1*hi+C2).*exp(C3.*hi.*hi+C4*hi);
    end

elseif pad == 2,

    logE=log10(Energy);
    if Energy >= 1000,
        lam=(Param_i(1,end)*hi+Param_i(2,end)).*exp(Param_i(3,end)*hi.*hi+Param_i(4,end)*hi);
    else
        C1=interp1(logE_i,Param_i(1,:),logE);
        C2=interp1(logE_i,Param_i(2,:),logE);
        C3=interp1(logE_i,Param_i(3,:),logE);
        C4=interp1(logE_i,Param_i(4,:),logE);
        lam=(C1*hi+C2).*exp(C3.*hi.*hi+C4*hi);
    end

end

end % energy_dissipation_Sergienko1993
