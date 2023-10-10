function [neStd,neCov] = ElSpec_error_estimate_ne2( ne0 , P , ne0Pcov , dt , A , photoprod , Ec , dE , alpha , integ_type )
%
% estimate standard deviation of the fitted electron density
% profile using the combined covariance matrix of ne0 and P. The
% function ElSpec_error_estimate_ne does the same error estimation,
% but uses the covariances of ne0 and Ie, but does not use their
% mutual covariances.
%
% neStd = ElSpec_error_estimate_ne( ne0 , ne00Cov , Ie , IeCov , dt
% , A , photoprod , Ec , dE , alpha , integ_type )
%
% IV 2017, 2018
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later


% density at the fit convergence point
nefit =  integrate_continuity_equation( ne0 , model_spectrum( P , Ec(:) ) , dt , A , photoprod , dE , alpha , integ_type );


% number of height gates
nh = length(ne0);

if any(any(isnan(ne0Pcov)))
    neStd = NaN(nh,1);
    neCov = NaN(nh,nh);
    return
end

% number of polynomial coefficients
nP = length(P);

% finite differences in P
dP = 1e-5*5.^[-(1:nP)];%1e-5*Ie + 1e3;

% finite differences in ne0
dne0 = ones(nh,1)*1e8;

% total number of "observations"
nm = nP + nh;

% the combined measurements
m = [ P' ; ne0 ];


% combined steps
dm  = [dP';dne0];

% linearize integrate_continuity_equation in vicinity of nefit
B = zeros(nh,nm);
for im = 1:nm
    mdiff = m;
    mdiff(im) = mdiff(im) + dm(im);
    B(:,im) = ( integrate_continuity_equation( mdiff(nP+1:end), ...
                                              model_spectrum(mdiff(1:nP),Ec(:)),dt,A,photoprod,dE, ...
                                              alpha,integ_type) - ...
                nefit ) / dm(im);
end

neCov = B * (ne0Pcov * B');

neStd = sqrt(diag(neCov));

end


