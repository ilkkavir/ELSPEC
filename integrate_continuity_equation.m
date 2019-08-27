function ne = integrate_continuity_equation(ne0,s,dt,A,dE,alpha,type)
%
% Integrate the electron continuity equation.
% The electron density, for cases where the convection can be
% ignored, can be calculated with a simplified continuity equation:
% 
% dne(z,t)/dt = -alpha(z,t)*ne(z,t)^2 + q(z,t)
% 
% if the recombination is dominated by dissociate recombination,
% i.e. recombination with molecular ions.
% 
% Usage: ne = integrate_continuity_equation(ne0,s,dt,A,alpha)
%
% INPUT:
%   ne0   electron density profile [m^-3] at beginning of the time step
%   s     electron precipitation spectrum [m^-2s eV^-1]
%   dt    time step [s]
%   A     a matrix of ion production profiles [m^-1]
%   dE    widths of the energy bins [eV]
%   alpha effective recombination coefficients [m^-3 s^-1]
%   type  a string telling which approximation to use:
%         'endNe'        do not integrate, but return ne after dt
%         'integrate'    actual integral value
%         'equilibrium'  solve ne from A*s = alpha*ne^2
%           (i.e. assume photochemical equilibrium, dne/dt = 0)
%         'linearend'    assume a simple linear model, production
%                        rate q and loss rate alpha*ne0^2, return
%                        ne after dt (do not integrate).
%         'linearint'    assume a simple linear model, production
%                        rate q and loss rate alpha*ne0^2,
%                        integrate ne over dt
%
%   ne0, A and alpha must be for the same altitudes.
%   s and A must be for the same energies
%
% OUTPUT:
%   ne    integral of electron density profile [m^-3] over time step dt
%
% I. Virtanen 2016, 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

% ion production rate
q = A*(s.*dE);

switch lower(type)
  case 'endne'
    ne = ( ne0 + sqrt(q./alpha).*tanh(sqrt(alpha.*q).*dt) ) ...
         ./ ( 1 + sqrt(alpha./q).*ne0.* tanh(sqrt(alpha.*q).*dt));
  case 'integrate'
    ne = ((log( sqrt(q) + sqrt(alpha).*ne0.*tanh(sqrt(alpha.*q).*dt)) - ...
           log(tanh(sqrt(alpha.*q).*dt) + 1) + sqrt(alpha.*q).*dt)./alpha ...
          - log(sqrt(q))./alpha ) / dt;
    izero = q <= 1e-100;
    if any(izero)
        ne(izero) = log( alpha(izero) .* ne0(izero) .* dt + 1) ./ ...
            alpha(izero) / dt;
    end
  case 'equilibrium'
    ne = sqrt(q./alpha);
  case 'linearend'
    % a simple linear model, production and loss rates are constant
    % during dt, return the end-point value...
    ne = ne0 + ( q - alpha.*ne0.^2) * dt;
    ne(ne<1) = 1;
  case 'linearint'
    % a simple linear model, production and loss rates are constant
    % during dt, return average over dt
    ne = ne0 + .5*( q - alpha.*ne0.^2) * dt;
    ne(ne<1) = 1;
  otherwise
    error(['Error. \nUnknown type %s in ' ...
           'integrate_continuity_equation.\nAccepted values ' ...
           'are ''','endNe''',', ''','integrate''',' and ''','equilibrium''','.'],type)

end

end

%% analytic solution from matlab
%%
%% notice that the form returned by pretty(ne_t) is numerically
%% unstable and a different form is used in the code!
%%
%%>> syms  q alpha ne0 ne(t) nei(t)
%%>> ne_t = dsolve(diff(ne)==q-alpha*ne^2,ne(0) == ne0);
%%>> nei_t = dsolve(diff(nei)==ne_t);
%%>> pretty(ne_t)
%%            /                     /          / sqrt(alpha) ne0 \\\
%%            |                     |     atanh| --------------- |||
%%            |                     |          \     sqrt(q)     /||
%%sqrt(q) tanh| sqrt(alpha) sqrt(q) | t + ------------------------||
%%            \                     \        sqrt(alpha) sqrt(q)  //
%%--------------------------------------------------------------------
%%                             sqrt(alpha)
%%
%%>> pretty(nei_t)
%%C5 + (log(sqrt(q) + sqrt(alpha) ne0
%%
%%   tanh(sqrt(alpha) sqrt(q) t)) - log(tanh(sqrt(alpha) sqrt(q)
%%
%%  t) + 1) + sqrt(alpha) sqrt(q) t)/alpha
%%
%%>> pretty(subs(nei_t,'t',0))
%%     log(sqrt(q))
%%C5 + ------------
%%         alpha
%%
