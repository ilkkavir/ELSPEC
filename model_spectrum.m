function s = model_spectrum(X0,E)
% energy spectrum of differential electron flux
% Model spectrum of the form E.*exp(X0(1)+X0(2)*E+X0(3)*E.^2,...)
%
% s = model_spectrum(X0,E)
%
% INPUT:
%  X0  coefficients of the polynomial
%  E   energy grid points
%
% OUTPUT:
%  s   the model spectrum
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

    persistent Epowers
    if isempty(Epowers) | size(Epowers,2) < numel(X0) | size(Epowers,1)~=numel(E)
        Epowers = ones(numel(E),numel(X0));
        for k = 1:numel(X0)
            Epowers(:,k) = (E(:)/1e4).^(k-1);
        end
    end

    
    % create the polynomial
    p = E(:).*0;
    for k=1:length(X0)
        %p = p + X0(k).*(E(:)/1e4).^(k-1);
        p = p + X0(k).*Epowers(:,k);
    end
    s = E(:).*exp( p(:) );

    
end
