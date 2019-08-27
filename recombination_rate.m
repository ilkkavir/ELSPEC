function [alpha1,alpha2,alpha3] = recombination_rate(Te,type)
%
% Recombination rates as function of electron temperature
% for O2+, N2+, and NO+.
%
% [alpha1,alpha2,alpha3] = recombination_rate(Te)
%
% INPUT:
%   Te     Electron temperature(s) [K]
%   type   a string defining which coefficients to use:
%          'Rees'      Rees (1989)
%          'SheehanGr' Sheehan and St.-Maurice (2004), ground state
%          'SheehanEx' Sheehan and St.-Maurice (2004),
%                      vibrationally excited states
%
% OUTPUT
%   alpha1 recombination rate for O2+ [m^3s^-1]
%   alpha2 recombination rate for N2+ [m^3s^-1]
%   alpha3 recombination rate for NO+ [m^3s^-1]
%
% IV 2016
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later



switch lower(type)
  case 'rees'

    %  alpha1: O2+ + e- -> O(1D) + O(eP); 5Ev
    alpha1 = 1.9e-13 * (Te/300).^(-.5);

    % alpha2: N2+ + e- -> N(2D) + N(2D); 1.04eV
    alpha2 = 1.8e-13 * (Te/300).^(-.39);

    % alpha3: NO+ + e- -> O + N(4S) + N(2D); 2.75eV
    alpha3 = 4.2e-13 * (Te/300).^(-.85);

  case 'sheehangr'

    %  alpha1: O2+ + e- -> O + O
    alpha1 = 1.95e-13 * (Te/300).^(-.7);

    % alpha2: N2+ + e- -> N + N
    alpha2 = 2.2e-13 * (Te/300).^(-.39);

    % alpha3: NO+ + e- -> O + N
    alpha3 = 3.5e-13 * (Te/300).^(-.69);


  case 'sheehanex'

    %  alpha1: O2+ + e- -> O + O
    alpha1 = 0.9e-13 * (Te/300).^(-.49);

    % alpha2: N2+ + e- -> N + N
    alpha2 = 1.5e-13 * (Te/300).^(-.39);

    % alpha3: NO+ + e- -> O + N
    alpha3 = 1.0e-13 * (Te/300).^(-.48);

  otherwise
    error(['Error. Unknown type %s in ' ...
           'recombination_rate.\nAccepted values ' ...
           'are ''','Rees''',', ''','SheehanGr''',', and ''','SheehanEx''','.'],type)



end

end