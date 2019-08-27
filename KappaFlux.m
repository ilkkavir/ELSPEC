function S = KappaFlux(E0,Q0,kappa,E)
%
% Differential number flux from kappa distribution
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

S = Q0/(2*E0.^3)*(kappa-1)*(kappa-2)/kappa^2*E.*(1+E/(kappa*E0)).^(-kappa-1);

end