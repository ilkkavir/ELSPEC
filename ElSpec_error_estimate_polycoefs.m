function  covar = ElSpec_error_estimate_polycoefs( pp , ppstd , ne0 , ne0cov, A, alpha , dt , Ec , dE , integtype , IePrior , stdPrior , polycoefs,nmeaseff,eType,ewidth)
%
% Estimate posterior covariance of the polynomial coefficients and ne0.
% Consider both ne0 and polycoefs as unknowns, and calculate
% Gaussian error estimates from finite-difference approximation
% of the Hessian matrix.
%
% covar = ElSpec_error_estimate_polycoefs( pp , ppstd , ne0 , ne0cov, A, alpha , dt , Ec , dE , integ_type_tt , IePrior , stdPrior , polycoefs , nmeaseff , etype , ewidth)
%
%
% IV 2017
% GB 2022
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later


% number of unknowns in the original fit
k = length(polycoefs);

% number of measurements in the original fit
n = length(pp(:));

% number of unknowns in error estimation
nx = k + length(ne0);

% hessian matrix
H = zeros(nx,nx);

% the parameter vector
param = [polycoefs , ne0'];

% the finite differences
dpar = [1e-5*5.^[-(1:k)],1e7*ones(1,nx-k)];

% inverse of the ne0 covariance
inescale = sqrt(diag(ne0cov));
ne0covs = ne0cov;
for nn=1:length(ne0)
    ne0covs(nn,:) = ne0covs(nn,:) / inescale(nn);
    ne0covs(:,nn) = ne0covs(:,nn) / inescale(nn);
    % the densities are strongly correlated, need to regularize to
    % actually get an inverse...
    ne0covs(nn,nn) = ne0covs(nn,nn) + 1e-10;
end
%ine0cov = inv(ne0covs);
ine0cov = pinv(ne0covs);
for nn=1:length(ne0)
    ine0cov(nn,:) = ine0cov(nn,:) / inescale(nn);
    ine0cov(:,nn) = ine0cov(:,nn) / inescale(nn);
end

% upper triangular part of the Hessian matrix (or, actually, -2*Hessian..), without the main diagonal
% the sum-of-squares is calculated by means of subtracting the AICc
% contribution from output of ElSpec_fitfun. This way we can use
% the exact same routines as in the iteration, and thus avoid
% bugs...
% the function does not account for error in ne0, but its effect
% will cancel out in the cross-terms anyway
for nn = 1:k%nx-1
    for ll = nn+1:nx

        % chi-squared at points par([nn,ll]) + dpar([nn,ll])
        par = param;
        par([nn,ll]) = par([nn,ll]) + dpar([nn,ll]);
        H(nn,ll) = H(nn,ll) + ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

        % chi-squared at points par([nn,ll]) - dpar([nn,ll])
        par = param;
        par([nn,ll]) = par([nn,ll]) - dpar([nn,ll]);
        H(nn,ll) = H(nn,ll) + ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

        %
        par = param;
        par([nn,ll]) = par([nn,ll]) + [dpar(nn) , -dpar(ll)];
        H(nn,ll) = H(nn,ll) - ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

        %
        par = param;
        par([nn,ll]) = par([nn,ll]) + [-dpar(nn) , dpar(ll)];
        H(nn,ll) = H(nn,ll) - ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

        % normalization
        H(nn,ll) = H(nn,ll) / (4 * dpar(nn) * dpar(ll));

    end
end

% the upper triangular part
H = H + H';

% chi-squared at the iteration end point
par = param;
chisqr0 = ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

% the main diagonal
for nn=1:k%nx
    par = param;
    par(nn) = par(nn) + dpar(nn);
    H(nn,nn) = H(nn,nn) + ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

    par = param;
    par(nn) = par(nn) - dpar(nn);
    H(nn,nn) = H(nn,nn) + ElSpec_fitfun(par(1:k),pp,ppstd,par(k+1:end),A,alpha,dt,Ec,dE,integtype,IePrior,stdPrior,nmeaseff,eType,ewidth);

    H(nn,nn) = H(nn,nn) - 2*chisqr0;

    H(nn,nn) = H(nn,nn) / dpar(nn)^2;

    % try to survive from numerical problems
    if H(nn,nn) <= 0
        H(nn,nn) = 1;
        %        disp('negitave element in diagonal of a Hessian, replaced with 1')
    end

end


% add contribution from the electron density covariance
H(k+1:end,k+1:end) = 2*ine0cov;

% normalize to avoid problems in inversion
H2 = .5*H;
hscales = sqrt(diag(H2));
for nn=1:nx
    H2(nn,:) = H2(nn,:)/hscales(nn);
    H2(:,nn) = H2(:,nn)/hscales(nn);
end

if any(any(isinf(H2))) | any(any(isnan(H2)))
    covar = H2 .* NaN;
    %    disp('Error estimation failed...')
    return
end

%covar = inv(H2);
covar = pinv(H2);

ntry = 0;
regconst = 1e-10;
while any(diag(covar)<0) | any(any(isnan(covar))) | any(any(imag(covar)~=0))
    for nn=1:nx
        H2(nn,nn) = H2(nn,nn) + regconst*10^ntry;
    end
    %    covar = inv(H2);
    try
        covar = pinv(H2);
    catch
        covar = H2.*NaN;
        return
    end
    ntry = ntry + 1;
    if ntry > 1000
        break
    end
end

%if ntry > 0
    %    disp('Regularization applied when inverting the Hessian matrix...')
    %    disp(ntry)
    %    disp(any(any(isnan(H2))))
    %    disp(any(diag(covar)<0))
    %    disp(any(any(isnan(covar))))
    %    disp(any(any(imag(covar)~=0)))
    %end

for nn=1:nx
    covar(nn,:) = covar(nn,:)/hscales(nn);
    covar(:,nn) = covar(:,nn)/hscales(nn);
end

% if any(sqrt(diag(covar))<1e-10)
%     disp('small variance!!')
%     disp(param(1:k))
%     disp(dpar(1:k))
% end

%imagesc(covar)
%drawnow
%pause(.5)
%disp('Elspec_error_estimate_polycoefs')
%disp(covar)
end

