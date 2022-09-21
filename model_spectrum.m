function s = model_spectrum(X0,E,S_type)
% energy spectrum of differential electron flux
% Model spectrum of the form E.*exp(X0(1)+X0(2)*E+X0(3)*E.^2,...)
%
% s = model_spectrum(X0,E)
%
% INPUT:
%  X0     coefficients of the polynomial
%  E      energy grid points
%  Stype  spectrum type empty or 'p' for the polynomial model, 'ger' for ??
%
% OUTPUT:
%  s   the model spectrum
%
% 
% IV 2017
% BG 2022
% 
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% create the polynomial
  if nargin < 3 || strcmp(lower(S_type),'p')
    p = E.*0;
    for k=1:length(X0)
      p = p + X0(k).*(E/1e4).^(k-1);
    end
    
    s = E.*exp( p );
  elseif strcmp(lower(S_type),'ger')
    X = [1 2 0 1 1];
    X(1:numel(X0)) = X0;
    s = E.^X(4).*exp(X(1)-abs((E-X(3)*1e3)/(X(2)*1e3)).^X(5));
  end
end