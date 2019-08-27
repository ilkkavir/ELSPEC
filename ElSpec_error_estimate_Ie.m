function covIe = ElSpec_error_estimate_Ie( covP , P , Ec )
%
% Covariance matrix of flux estimates on energy grid points,
% based on  covariance estimate for the polynomial coefficients
%
% covIe = Elspec_error_estimate_Ie( covP , P , Ec )
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% the flux estimate at the fit convergence point
Ie = model_spectrum( P , Ec );

% number of energy bins
nE = length(Ec);

% any NaN's will produce NaN output anyway...
if any(any(isnan(covP)))
    covIe = NaN(nE,nE);
    return
end


% number of polynomial coefficients
k = length(P);

% finite differences of P
dP = 1e-5*5.^[-(1:k)];
%dP = ones(1,k)*1e-6;


% linearize model_spectrum in vicinity of Ie
A = zeros(nE,k);
for iP = 1:k
    Pdiff = P;
    Pdiff(iP) = Pdiff(iP) + dP(iP);
    A(:,iP) = (model_spectrum( Pdiff , Ec ) - Ie)/dP(iP);
end


covIe = A * ( covP * A' );

if any(any(imag(covIe))~=0)
    covIe = real(covIe);
    %    disp('ElSpec_error_estimate_Ie: removing imaginary part...')
end

if any(diag(covIe)<0)
    covIe = covIe - diag(ones(1,nE))*min(diag(covIe)-1e-7);
    %    disp('ElSpec_error_estimate_Ie: regularizing...')
end


%if any(diag(covIe)<0) | any(any(imag(covIe)~=0))
%    disp('ElSpec_error_estimate_Ie failed...')
%%    disp('covIe')
%disp(diag(covIe))
%disp(diag(A*A'))
%    disp(any(diag(covIe)<0))
%    disp(any(any(imag(covIe)~=0)))
%    disp(covP)
%    covIe = NaN(nE,nE);
%    return
%end

end