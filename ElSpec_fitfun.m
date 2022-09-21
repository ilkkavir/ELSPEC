function AIC = ElSpec_fitfun(X0,ne,stdne,ne0,A,alpha,dt,E,dE,integtype,IePrior,stdPrior,nmeas_eff,eType,eWidth,varargin)
%                             1  2     3   4 5     6  7 8  9        10      11       12        13    14     15       16    17
% The function to be minimized in the electron energy spectrum
% fit. Re-arranges the paremeter vector p in a more understanble
% form and calls the function update_AICc.
%
% AIC = ElSpec_fitfun(X0,ne,stdne,ne0,A,alpha,dt,E,integtype)
%
% INPUT:
%  X0         Coefficients of the polynomial
%  ne         measured electron density profile
%  stdne      standard deviations of ne
%  ne0        modeled electron density profile from the previous time step
%  A          Ion production profile matrix
%  alpha      effective recombination rates
%  dt         time step [s]
%  E          energies of the dense "model grid"
%  dE         widths of the energy bins
%  integtype  type input for integrate_continuity_equation
%             ('endNe','integrate',or 'equilibrium')
%  IePrior    Apriori flux
%  stdPrior   standard deviation of the apriori
%  eType      type of error distribution: 's', 'n', 't', or 'l'.
%  eWidht     width of transition from normal to long-tailed distribution in units of sigma
%
%
% OUTPUT:
%  AIC        corrected Akaike information criterion value
%
%
%
% IV 2017
% BG 2022
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later



% calculate the actual AICc value
if nargin > 15 && ~isempty(varargin{1})
  Ie_type = varargin{1};
  AIC = update_AICc(ne,stdne,ne0,A,alpha,dt,X0,E,dE,integtype,IePrior,stdPrior,nmeas_eff,eType,eWidth,Ie_type);
else
  AIC = update_AICc(ne,stdne,ne0,A,alpha,dt,X0,E,dE,integtype,IePrior,stdPrior,nmeas_eff,eType,eWidth);
end

%
% If we need this regularization, it should be move to update_AICc
%
% if size(X0,1) > 1
%   PenParVar = sum(diff(X0).^2,'all','omitnan');
%   AIC = AIC + PenParVar;
% end

end