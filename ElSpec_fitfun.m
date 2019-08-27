function AIC = ElSpec_fitfun(X0,ne,stdne,ne0,A,alpha,dt,E,dE,integtype,IePrior,stdPrior,nmeas_eff)
%
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
%
%
% OUTPUT:
%  AIC        corrected Akaike information criterion value
%
%
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later



% calculate the actual AICc value
AIC = update_AICc(ne,stdne,ne0,A,alpha,dt,X0,E,dE,integtype,IePrior,stdPrior,nmeas_eff);
%if isnan(AIC)
%    disp('AIC=NaN')
%end

end