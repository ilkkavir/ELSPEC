function [q,Ec,dE] = ion_production(E,z,nN2,nO2,nO,nAr,Tn,model)
%
% Ionozation profiles of monoenergetic electron precipitation
% using selected model.
%
%   [q,Ec,dE] = ion_production(E,z,nN2,nO2,nO,Tn)
%
%
%
% INPUT:
%   E     electron energies [eV] (1 x (nE+1) array)
%   z     altitudes [m] (nz x 1 array)
%   nN2   molecular Nitrogen density [m^-3] (nz x 1 array)
%   nO2   molecular Oxygen density [m^-3] (nz x 1 array)
%   nO    atomic oxygen density [m^-3] (nz x 1 array)
%   nAr   Argon density [m^-3] (nz x 1 array)
%   Tn    neutral temperature [K] (nz x 1 array)
%   model The model to use, either 'Fang' or 'Sergienko'.
%
% OUTPUT:
%   q   ionization profiles [m^-1] (nz x nE array)
%       at centres of the energy bins, i.e. energies (E(1)+E(2))/2,
%       (E(2)+E(3))/2 etc.
%   Ec  centre points of the energy bins [eV]
%   dE  widths of the energy bins [eV]
%
% Details:
%
%  The 'Fang' model uses the Fang et al. (2010) model in function
%  ion_production_Fang2010. See
%  help ion_production_Fang2010
%   for details.
%
%  The 'Sergienko' model uses the Sergienko et al. (1993) model in
%  function ion_production_Sergienko1993. See
%  help ion_production_Sergienko1993
%  for details.
%
%
% I. Virtanen 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% use persistent variables to avoid unnecessary re-calculation
persistent prev_input prev_q prev_Ec prev_dE

this_input = [E(:);z(:);nN2(:);nO2(:);nO(:);Tn(:);double(model)'];

% if the previous input was stored
if 0% ~isempty(prev_input)
    % the lengths may not match if 'model' was different
    if length(prev_input) == length(this_input)
        % if all inputs were identical, we already have the results
        if all(prev_input == this_input)
            q = prev_q;
            Ec = prev_Ec;
            dE = prev_dE;
            return
        end
    end
end

% call the actual model functions

switch lower(model(1:4))
  case 'fang'
    [q,Ec,dE] = ion_production_Fang2010(E,z,nN2,nO2,nO,nAr,Tn);
  case 'serg'
    [q,Ec,dE] = ion_production_Sergienko1993(E,z,nN2,nO2,nO,Tn);
  case 'rees'
    [q,Ec,dE] = ion_production_Rees1963(E,z,nN2,nO2,nO,Tn);
  otherwise
   error(['Unknowon ionization model "' , model , '". Use "Fang", ' ...
                       ' "Sergienko", or "Rees".'])
end

% make persistent copies of the input and output
prev_input = this_input;
prev_q = q;
prev_Ec = Ec;
prev_dE = dE;

end