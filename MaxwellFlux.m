function S = MaxwellFlux(E0,Q0,E)
%
% differential number flux from maxwellian energy distribution
%
% CALLING:
%   S = MaxwellFlux(E0,Q0,E)
% 
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

S = Q0/(2*E0.^3)*E.*exp(-E/E0);

end