function AIC = update_AICc(ne,stdne,ne0,A,photoprod,alpha,dt,X0,E,dE, ...
                           integtype,IePrior,stdPrior , nmeas)
%
% Calculate corrected Akaike information criterion value for model
% value of the  next electron density profile.
%
% ne = update_AICc(ne,stdne,ne0,A,alpha,dt,E0,Y0,gamma,E,integtype)
%
% INPUT:
%  ne         measured electron density profile
%  stdne      standard deviations of ne
%  ne0        modeled electron density profile from the previous time step
%  A          Ion production profile matrix
%  photoprod  Ion production by photoionization
%  alpha      effective recombination rates
%  dt         time step [s]
%  X0         coefficients of the polynomial
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
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% the high-resolution spectrum
s = model_spectrum( X0 , E );

% the updated model ne profile
if strcmp(integtype,'equilibrium')
    nemod = integrate_continuity_equation(ne0(:),s(:),dt(1),A,photoprod,dE(: ...
                                                      ),alpha(:),integtype);
    nemeas = ne(:,1);
    stdmeas = stdne(:,1);
else
    switch lower(integtype)
      case 'endne'
        integ_type_end = 'endNe';
      case 'integrate'
        integ_type_end = 'endNe';
      case 'linearend'
        integ_type_end = 'linearend';
      case 'linearint'
        integ_type_end = 'linearend';
      otherwise
        error(['Error. \nUnknown integtype %s in update_AICc' ] , type )
    end

    ninteg = length(dt);
    nemod = NaN(length(ne0),ninteg);
    neEnd = NaN(length(ne0),ninteg);
    nemod(:,1) = integrate_continuity_equation(ne0(:),s(:),dt(1),A,photoprod,dE(:),alpha(:),integtype);
    neEnd(:,1) = integrate_continuity_equation(ne0(:),s(:),dt(1),A,photoprod,dE(:),alpha(:),integ_type_end);
    if ninteg > 1
        for tt=2:ninteg
            nemod(:,tt) = integrate_continuity_equation(neEnd(:,tt-1),s(:),dt(tt),A,photoprod,dE(:),alpha(:),integtype);
            neEnd(:,tt) = integrate_continuity_equation(neEnd(:,tt-1),s(:),dt(tt),A,photoprod,dE(:),alpha(:),integ_type_end);
        end
    end
    nemeas = ne;
    stdmeas = stdne;
end

% the information criterion value
%AIC = AICc( nemeas(:) , stdmeas(:).^2 , nemod(:) , length(X0) , length(nemeas(:)) );
AIC = AICc( nemeas(:) , stdmeas(:).^2 , nemod(:) , length(X0) , nmeas );

% slight regularization for all coefficients
AIC = AIC + sum(X0.^2./1e5.*10.^[1:length(X0)]);

%% the spectrum should go to zero at the high energy end
%AIC = AIC + s(end)^2;

%% reject negative densities, these are possible only with the
%% linear model...
%AIC = AIC + any(nemod(:)<1)*1e10;


% try to damp oscillations by means of smoothing...
% allow larger variations when difference between previous
% solution and the new measured Ne is large..
%if ~isempty(Ieprev)
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)./1e18);
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)/mean(9e18*(ne0-ne).^2./stdne.^2));
%    %    AIC = AIC + sum(((Ieprev(:)-s(:)).^2)/(1e15*mean(abs(ne0.^4-ne.^4)./stdne.^4)));
%
%end

% try with a prior model (from an equilibrium solution...)
if ~isempty(IePrior)
    AIC = AIC + mean( ( s(:) - IePrior(:) ).^2 ./(stdPrior(:)).^2);
end

% tell the solver that flux at 100 eV is probably not huge...
%s100 = model_spectrum(X0,100);
%AIC = AIC + s100.^2/1e18;


end
