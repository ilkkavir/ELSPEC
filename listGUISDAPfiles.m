function [fpaths] = listGUISDAPfiles( ddir )
% Return list of .mat-files profduced by GUISDAP
% 
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later


% list mat files
ff = dir( fullfile( ddir , '*.mat') );
% GUISDAP output files have 8 digits and '.mat' in file name.
% The digits are seconds since beginning of year, so they are
% necessarily between 0 and 366*24*60*60=31622400
for k=length(ff):-1:1 % we will be removing elements, must go backwards
    if length(ff(k).name)~=12
        ff(k) = [];
    % str2num returns an empty matrix if its argument does not
    % represent a number
    elseif isempty( str2num(ff(k).name(1:8)) )
        ff(k) = [];
    % the number is seconds from beginning of the year and must
    % thus be larger than zero and smaller than 366*24*60*60
    elseif ( str2num(ff(k).name(1:8)) < 0 ) | ( str2num(ff(k).name(1:8)) > 31622400 )
        ff(k) = [];
    end
end

fpaths=cell(0);
% if nothing was left
if isempty(ff)
    return;
end

% full paths to data files
for k=1:length(ff)
    fpaths{k} = fullfile( ddir , ff(k).name );
end
