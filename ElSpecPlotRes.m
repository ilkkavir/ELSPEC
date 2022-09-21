function fignum = ElSpecPlotRes(ElSpecOut,varargin)
% ElSpecPlotRes - Plot ElSpec Output, plots electron-density residuals only
%
% fignum = ElSpecPlotRes( ElSpecOut , <fignum> , ... )
%
% INPUT:
%  ElSpecOut   an output list from ElSpec
%  fignum      optional number of the figure. A new figure is
%              created by default.
%  ...         optional parameters as name-value pairs:
%    ielim     differential number flux color scale limits, (log((eVm^2s)^-1)),
%              default [6 10]
%    ieelim    differential energy flux color scale limits
%              (log(eV/(eVm^2s))), default [10 14]
%    elim      energy axis limits for the differential flux plot
%              (keV), default [1 100]
%    reslim    Residual limits, default [-3 3]
%    plim      net energy flux limits (W/^2), default [0 50]
%    btime     plot start time as vector [yyyy mm dd HH MM SS],
%              default: first data point
%    etime     plot end time as vector [yyyy mm dd HH MM SS],
%              default: last data point
%    fluxtype  flux plot type, 'number', 'energy', or 'both'. If
%              'both' is selected, the net energy flux plot is
%              omitted. Default: 'energy'
%    emin      minimum energy to include in the FAC and net energy
%              flux integrations (eV). Default 1000. NOTICE: error
%              estimates cannot be calculated if 'emin' is changed
%              and ElSpecOut was calculated with
%              'saveiecov','false' in ElSpec.
%    cutgaps   logical, should white space be plotted on obvious
%              data gaps. Default true
%
% Modified by B Gustavsson 2022, from ElSpecPlot.m by I Virtanen:
% IV 2016, 2017, 2018, 2019
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


p = inputParser;

defaultFignum = -1;

defaultHlim = [ElSpecOut.h(1) ElSpecOut.h(end)];
checkHlim = @(x) (isnumeric(x) & (length(x)==2));

defaultReslim = [-3 3];
checkReslim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecOut.te(1) + ElSpecOut.ts(1))/2,'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecOut.te(end) + ElSpecOut.ts(end))/2,'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultCutgaps = 1;
checkCutgaps = @(x) (islogical(x)|isnumeric(x));

addRequired(p,'ElSpecOut',@isstruct);
addOptional(p,'fignum',defaultFignum,@isnumeric);
addParameter(p,'hlim',defaultHlim,checkHlim)
addParameter(p,'reslim',defaultReslim,checkReslim)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'cutgaps',defaultCutgaps,checkCutgaps)
parse(p,ElSpecOut,varargin{:})


% put white space on data gaps
if p.Results.cutgaps
    mediandt = median(diff(ElSpecOut.ts));
    gaps = diff(ElSpecOut.ts) > (10*mediandt);
    rminds = [find(gaps) ; (find(gaps)+1)];
    ElSpecOut.pp(:,rminds) = NaN;
    ElSpecOut.ne(:,rminds) = NaN;
    ElSpecOut.Ie(:,rminds) = NaN;
    ElSpecOut.FAC(rminds) = NaN;
    ElSpecOut.FACstd(rminds) = NaN;
    ElSpecOut.Pe(rminds) = NaN;
    ElSpecOut.PeStd(rminds) = NaN;
    ElSpecOut.chisqr(rminds) = NaN;
end


% centres of time-bins as matlab datenum
tt = datenum(datetime((ElSpecOut.te + ElSpecOut.ts)/2,'ConvertFrom','posixtime'));
% start points of time-bins for the pcolor plots
ts = datenum(datetime(ElSpecOut.ts,'ConvertFrom','posixtime'));
% end points of time-bins for the pcolor plots
te = datenum(datetime(ElSpecOut.te,'ConvertFrom','posixtime'));
% heighs
hh = ElSpecOut.h;

% Electron densities
ne = ElSpecOut.ne;
% power profiles with negatives removed
pp = ElSpecOut.pp;
pp(pp<1) = 1;
ppstd = ElSpecOut.ppstd;

fignum = p.Results.fignum;

if fignum>0
    figure(fignum);
    figpos = get(gcf,'Position');
else
    figtmp = figure;
    fignum = figtmp.Number;
    figpos = [0 0 21/29.7*690 690];
end
figpos(3:4) = [21/29.7*690 690*.6];
set(gcf,'Position',figpos,'PaperPositionMode','Auto');

h4 = gca;

pcolor([ts(:)',te(end)],hh,(pp(:,[1:end,end]) - ne(:,[1:end,end]))./ppstd(:,[1:end,end])),shading flat
caxis(p.Results.reslim)
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
ylabel('Height [km]')
cbh4 = colorbar;
ylabel(cbh4,{'\Delta N_e/\sigma_{n_e}'})

datetick(gca,'x',13,'keeplimits','keepticks')

cbsize1 = get(cbh4,'Position');

set(cbh4,'fontsize',12)

set(h4,'fontsize',12)

set(h4,'layer','top','colormap',redblue)

drawnow


end