function ss = AICc( m , v , y , k  , n )
%
% Corrected Akaike information criterion.
%
% Usage: ss = AICc(m,v,y,k,n)
%
% Akaike information criterion corrected for finite sample size.
% The correction assumes that the model is univariate, linear, and
% has normally-distributed residuals.
%
% INPUT:
%   m    vector of measurements
%   v    variances of the measurements
%   y    direct theory values
%   k    number of parameters
%   n    sample size
%
% OUTPUT:
%   ss   Corrected Akaike information criterion
%
%
% IV 2016
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later


ss = ( sum( ( m - y ).^2 ./ v ) ) + 2*k + 2*k*(k+1)/(n-k-1);

end