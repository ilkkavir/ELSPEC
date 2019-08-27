function [q,Ec,dE] = ion_production_Fang2010(E,z,nN2,nO2,nO,nAr,Tn)
%   Ionization profiles of monoenergetic electron precipitation
%   according to Fang et al [2010]. This function calculates the
%   ionization profiles for electron precipitation into the Earths
%   atmosphere using a parameterization that are accurate in the
%   energy range 100 eV - 1 MeV, both in total ionization rate
%   (error < +/-4%), peak ionization rate (error < 10%),altitude of
%   peak ionization (error < 2 km E>500 eV, < 5 km, E < 500 eV).
%
%   Based on the function q_e_z_Fang2 by B. Gustavsson.
%
%   [q,Ec,dE] = ion_production_Fang2010(E,z,nN2,nO2,nO,Tn)
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

% Altitude varying gravitational acceleration
g = 9.82*(Re./(Re+z')).^2;

% mass density (kg/m^3)
ro = (nO*8*(m_n+m_p+m_e) + nO2*2*8*(m_n+m_p+m_e) + nN2*2*7.5*(m_n+m_p+m_e)) ...
     + nAr*(18*m_p+22*m_n+18*m_e);

% Average molecular mass (kg)
m = ro./(nO + nO2 + nN2 + nAr);
q = zeros(nz,nE);
for iE = 1:nE,

    E0 = Ec(iE)/1e3;                         % Rescale electron energy to keV
    H = (kB*Tn./(m.*g))*100;                 % scale height [cm]
    y = 2./(E0).*((ro/1000).*H/6e-6).^(0.7); % "normalized atmospheric column mass"
    dEion = 0.035;                           % Average energy per ion-electron-pair [keV]
    % ionization profile in m^-1
    q(:,iE) = energy_dissipation_Fang2010(y,E0)/dEion./H*E0*100;

end


end %ion_production_Fang2010

function f = energy_dissipation_Fang2010(y,E)
% normalized eneryg dissipation function
%
% Calling:
%  f = energy_dissipation_Fang2010(y,E)
% Input:
%  y - normalized atmospheric column mass, double array [nZ x 1]
%  E - energy, double scalar (keV)
% Output:
%  f - normalized energy deposition

f = Ci_Fang2010(E,1).*y.^Ci_Fang2010(E,2).*exp(-Ci_Fang2010(E,3).*y.^Ci_Fang2010(E,4)) + ...
    Ci_Fang2010(E,5).*y.^Ci_Fang2010(E,6).*exp(-Ci_Fang2010(E,7).*y.^Ci_Fang2010(E,8));


end % energy_dissipation_Fang2010

function c_i = Ci_Fang2010(E,i)
% Ci_Fang2010 - exponential of cubic polynomial
%
% Calling:
%  c_i = Ci_Fang2010(E,i)
% Input:
%  E - energy, double scalar (keV)
%  i - polynom-index, integer scalar [1-8], the Fang-et.-al. paper
%      calculates the normalized energy deposition as an expresion
%      which contains eight energy-dependent factors that are
%      exponentials of cubic polynomials.
% Output:
%  c_i - exp(sum(P(i,:)*log(E).^[0:3]));

j = [0          1           2           3];
Pij = [...
    1.24616e+0  1.45903e+0 -2.42269e-1  5.95459e-2   % i = 1
    2.23976e+0 -4.22918e-7  1.36458e-2  2.53332e-3   % i = 2
    1.41754e+0  1.44597e-1  1.70433e-2  6.39717e-4   % i = 3
    2.48775e-1 -1.50890e-1  6.30894e-9  1.23707e-3   % i = 4
   -4.65119e-1 -1.05081e-1 -8.95701e-2  1.22450e-2   % i = 5
    3.86019e-1  1.75430e-3 -7.42960e-4  4.60881e-4   % i = 6
   -6.45454e-1  8.49555e-4 -4.28581e-2 -2.99302e-3   % i = 7
    9.48930e-1  1.97385e-1 -2.50660e-3 -2.06938e-3]; % i = 8


c_i = 0*E;
for j1 = j,
  c_i = c_i + Pij(i,j1+1)*log(E).^j1;
end
c_i = exp(c_i);

end %Ci_Fang2010