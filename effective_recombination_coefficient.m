function alpha = effective_recombination_coefficient(h,Te,O2p,N2p,NOp,type)
%
% Effective electron recombination coefficient, i.e. alpha in the equation
% dne/dt = q - alpha*ne^2.
% 
% alpha = effective_recombination_coefficient( h , Te , O2p , N2p , NOp , type )
%
% INPUT:
%   h    altitude vector [km]
%   Te   electron temperatures [K]
%   O2p  O2 ion fractions
%   N2p  N2 ion fractions
%   NOp  NO ion fractions
%   type string defining the type of coefficient calculation:
%        'Rees'     calculate alpha from Te and composition using
%                   formulae of Rees (1989) and the input ion
%                   composition.
%        'ReesO2+'   Rees (1989), assuming pure O2+
%        'ReesN2+'   Rees (1989), assuming pure N2+
%        'ReesNO+'   Rees (1989), assuming pure NO+
%        'SheehanGr' Sheehan and St.-Maurice (2004), ground state,
%                    input composition
%        'SheehanGrO2+'   Sheehan and St.-Maurice (2004) ground state, assuming pure O2+
%        'SheehanGrN2+'   Sheehan and St.-Maurice (2004) ground state, assuming pure N2+
%        'SheehanGrNO+'   Sheehan and St.-Maurice (2004) ground state, assuming pure NO+
%        'SheehanEx' Sheehan and St.-Maurice (2004), vibrationally
%                    excited, input composition
%        'delPozo1' alpha as function of height only, eq. 13 from
%                   del Pozo et al. (1997):
%                   alpha = 2.5e-12*exp(-.0242*h) +
%                   1.63e5*exp(-.524*h);
%        'delPozo2' alpha as function of height only, an even
%                   simpler equation from del Pozo et al.:
%                   alpha =  = 2.5e-12*exp(-.0195*h);
%
% OUTPUT:
%  alpha effective recombination coefficients [m^3 s^-1] at heights h.
%
%  I. Virtanen 2016
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later


% BG suggestions:
% use default type? And use that as fallback for invalid cases?
if nargin < 6 || isempty(type)
  type = 'rees';
end
% And use that as fallback for invalid cases?

switch lower(type)
  case 'rees'
    [a1,a2,a3] = recombination_rate( Te , 'Rees' );
    alpha = sum( [a1 a2 a3] .* [O2p N2p NOp] , 2 );
    % make sure that sum of ion abundances is 1
    alpha = alpha ./ sum([O2p N2p NOp] , 2 );
  case 'reeso2+'
    [a1,a2,a3] = recombination_rate( Te , 'Rees' );
    alpha = a1;
  case 'reesn2+'
    [a1,a2,a3] = recombination_rate( Te , 'Rees' );
    alpha = a2;
  case 'reesno+'
    [a1,a2,a3] = recombination_rate( Te , 'Rees' );
    alpha = a3;
  case 'delpozo1'
    alpha = 2.5e-12*exp(-.0242*h) + 1.63e5*exp(-.524*h);
  case 'delpozo2'
    alpha = 2.5e-12*exp(-.0195*h);
  case 'const'
    alpha = h.*0+4e-13;
  case 'sheehangr'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanGr' );
    alpha = sum( [a1 a2 a3] .* [O2p N2p NOp] , 2 );
    % make sure that sum of ion abundances is 1
    alpha = alpha ./ sum([O2p N2p NOp] , 2 );
  case 'sheehagro2+'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanGr' );
    alpha = a1;
  case 'sheehangrn2+'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanGr' );
    alpha = a2;
  case 'sheehangrno+'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanGr' );
    alpha = a3;
  case 'sheehangrflipchem'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanGr' );
    alpha = sum( [a1 a2 a3] .* [O2p N2p NOp] , 2 );
    % make sure that sum of ion abundances is 1
    alpha = alpha ./ sum([O2p N2p NOp] , 2 );    
  case 'sheehanex'
    [a1,a2,a3] = recombination_rate( Te , 'SheehanEx' );
    alpha = sum( [a1 a2 a3] .* [O2p N2p NOp] , 2 );
    % make sure that sum of ion abundances is 1
    alpha = alpha ./ sum([O2p N2p NOp] , 2 );
  otherwise
   % And use that as fallback for invalid cases instead of erroring
   % out and change this to a warning, or raw disp (so no warning
   % off can silence the message)?
   
    error(['Error. Unknown type %s in ' ...
           'effective_recombination_coefficient.\nAccepted values ' ...
           'are ''','Rees''',', ''','ReesO2+''',', ''','ReesN2+''',', ' ...
           '''','ReesNO+''',', ''','delPozo1''',', ''','delPozo2''',', ''','SheehanGr''',', ''','SheehanEx''',', and ''','const''','.'],type)



end

end