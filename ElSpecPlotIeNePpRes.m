function fignum = ElSpecPlotIeNePpRes(ElSpecOut,varargin)
% ElSpecPlotIeNePpRes - Plot ElSpec Output, this function plots only the
% energy spectrum, the modeled electron-density, the observed
% electron-density and the normalised residuals
%
% fignum = ElSpecPlot( ElSpecOut , <fignum> , ... )
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
%
% IV 2016, 2017, 2018, 2019
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


p = inputParser;

defaultFignum = -1;

defaultHlim = [ElSpecOut.h(1) ElSpecOut.h(end)];
checkHlim = @(x) (isnumeric(x) & (length(x)==2));

defaultIelim = [6 10];
checkIelim = @(x) (isnumeric(x) & (length(x)==2));

defaultIEelim = [10 14];
checkIEelim = @(x) (isnumeric(x) & (length(x)==2));

defaultNelim = [10 12.25];
checkNelim = @(x) (isnumeric(x) & (length(x)==2));

defaultElim = [1 100];%[min(ElSpecOut.Ec) max(ElSpecOut.Ec)]./1000;
checkElim = @(x) (isnumeric(x) & (length(x)==2));

defaultReslim = [-3 3];
checkReslim = @(x) (isnumeric(x) & (length(x)==2));

defaultPlim = [0 50];
checkPlim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecOut.te(1) + ElSpecOut.ts(1))/2,'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecOut.te(end) + ElSpecOut.ts(end))/2,'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultFluxtype = 'energy';
validFluxtype = {'number','energy','both'};
checkFluxtype = @(x) any(validatestring(x,validFluxtype));

defaultNeplot = 'log';
validNeplot = {'log','linear'};
checkNeplot = @(x) any(validatestring(x,validNeplot));

defaultCutgaps = 1;
checkCutgaps = @(x) (islogical(x)|isnumeric(x));

if isfield(ElSpecOut,'Emin')
    defaultEmin = ElSpecOut.Emin;
else
    defaultEmin = 1e3;
end
checkEmin = @(x) (isnumeric(x) & length(x)==1 );

addRequired(p,'ElSpecOut',@isstruct);
addOptional(p,'fignum',defaultFignum,@isnumeric);
addParameter(p,'hlim',defaultHlim,checkHlim)
addParameter(p,'ielim',defaultIelim,checkIelim)
addParameter(p,'ieelim',defaultIEelim,checkIEelim)
addParameter(p,'elim',defaultElim,checkElim)
addParameter(p,'nelim',defaultNelim,checkNelim)
addParameter(p,'reslim',defaultReslim,checkReslim)
addParameter(p,'plim',defaultPlim,checkPlim)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'fluxtype',defaultFluxtype,checkFluxtype)
addParameter(p,'emin',defaultEmin,checkEmin)
addParameter(p,'cutgaps',defaultCutgaps,checkCutgaps)
addParameter(p,'neplot',defaultNeplot,checkNeplot)
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
% energies (keV)
EE = ElSpecOut.Ec/1000;

% Electron densities
ne = ElSpecOut.ne;
% power profiles with negatives removed
pp = ElSpecOut.pp;
pp(pp<1) = 1;
ppstd = ElSpecOut.ppstd;

% number vs energy flux
switch lower(p.Results.fluxtype)

  case 'number'
    Ie = ElSpecOut.Ie;
    Ielabel={'I_e [(eVm^2s)^{-1}]'};
    ielim = p.Results.ielim;
  case 'energy'
    Ie = ElSpecOut.Ie;
    for iiE = 1:length(ElSpecOut.Ec)
        Ie(iiE,:) = Ie(iiE,:).*ElSpecOut.Ec(iiE);
    end
    Ielabel={'I_e [eV/(eVm^2s)]'};
    ielim = p.Results.ieelim;
  case 'both'
    Ie = ElSpecOut.Ie;
    Ielabel={'I_e [(eVm^2s)^{-1}]'};
    IEe = ElSpecOut.Ie;
    for iiE = 1:length(ElSpecOut.Ec)
        IEe(iiE,:) = IEe(iiE,:).*ElSpecOut.Ec(iiE);
    end
    IEelabel={'I_e [eV/(eVm^2s)]'};
    ielim = p.Results.ielim;
  otherwise
    error(['Unknown flux type' , p.Results.fluxtype])
end


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
try
  colormap(turbo);
catch
  colormap(jet);
end

h1 = subplot(4,1,1);
pcolor([ts(:)',te(end)],EE,log10(Ie(:,[1:end,end]))),shading flat
set(h1,'yscale','log')
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
caxis(ielim)
ylabel('Energy [keV]')
cbh1 = colorbar;
ylabel(cbh1,Ielabel)

h2 = subplot(4,1,2);
if strcmp(p.Results.neplot,'log')
  pcolor([ts(:)',te(end)],hh,log10(max(0,ne(:,[1:end,end])))),shading flat
  caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
  pcolor([ts(:)',te(end)],hh,ElSpecOut.ne(:,[1:end,end])),shading flat
  caxis(10.^(p.Results.nelim))
else
  error(['Unknown neplot ',neplot])
end
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.hlim)
cbh2 = colorbar;
ylabel('Height [km]')
ylabel(cbh2,{'N_e [m^{-3}] (model)'})


grid on


h3 = subplot(4,1,3);
if strcmp(p.Results.neplot,'log')
  pcolor([ts(:)',te(end)],hh,log10(max(0,pp(:,[1:end,end])))),shading flat
  caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
  pcolor([ts(:)',te(end)],hh,pp(:,[1:end,end])),shading flat
  caxis(10.^(p.Results.nelim))
else
error(['Unknown neplot ',neplot])
end
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
ylabel('Height [km]')
cbh3 = colorbar;
ylabel(cbh3,{'N_e [m^{-3}] (meas.)'})


h4 = subplot(4,1,4);
pcolor([ts(:)',te(end)],hh,-(ne(:,[1:end,end]) - pp(:,[1:end,end]))./ppstd(:,[1:end,end])),shading flat
caxis(p.Results.reslim)
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
ylabel('Height [km]')
cbh4 = colorbar;
ylabel(cbh4,{'\Delta N_e/\sigma_{n_e}'})




drawnow % this is important to update everything before setting the sizes...

pos1 = get(h1,'Position');
pos2 = get(h2,'Position');
pos3 = get(h3,'Position');
pos4 = get(h4,'Position');
set(h1, 'Position', [pos1(1:2)-[0,.06] [.85 1.4].*pos3(3:4)]);
set(h2, 'Position', [pos2(1:2)-[0,.04] [.85 1.4].*pos3(3:4)]);
set(h3, 'Position', [pos3(1:2)-[0,.02] [.85 1.4].*pos3(3:4)]);
set(h4, 'Position', [pos4(1:2)-[0,.02] [.85 1.4].*pos3(3:4)]);

linkaxes([h1 h2 h3 h4],'x')

eticks=[.5 1 2 4 8 16 32 64 128 256 512 1024];
set(h1,'YTick',eticks(eticks<p.Results.elim(2)))
set(h1,'TickDir','both')
set(h2,'TickDir','both')

set(h3,'TickDir','both')
set(h4,'TickDir','both')

set(h2,'XTick',get(h1,'XTick'))
set(h3,'XTick',get(h1,'XTick'))
set(h4,'XTick',get(h1,'XTick'))
datetick(h2,'x',13,'keeplimits','keepticks')
datetick(h3,'x',13,'keeplimits','keepticks')
datetick(h4,'x',13,'keeplimits','keepticks')
set(h1,'XTickLabel','')
set(h2,'XTickLabel','')
set(h3,'XTickLabel','')
set(h1,'LineWidth',1)
set(h2,'LineWidth',1)
set(h3,'LineWidth',1)
set(h4,'LineWidth',1)

cbsize1=get(cbh1,'Position');
% $$$ set(cbh1,'Position',cbsize1.*[1,1.02,1,.8]);
set(cbh1,'fontsize',12)

set(h1,'fontsize',12)
set(h2,'fontsize',12)
set(h3,'fontsize',12)
set(h4,'fontsize',12)

set(cbh1,'YTick',5:13)
set(cbh1,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})
set(cbh2,'YTick',5:13)
set(cbh2,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})
set(cbh3,'YTick',5:13)
set(cbh3,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})


set(h1,'layer','top')
set(cbh1,'LineWidth',1,'TickDir','both')
set(h2,'layer','top')
set(h3,'layer','top')
set(h4,'layer','top','colormap',redblue)

tstr1 = datestr(datenum(datetime((ElSpecOut.te(1)),'ConvertFrom', ...
                                 'posixtime')),29);
tstr2 = datestr(datenum(datetime((ElSpecOut.te(end)),'ConvertFrom', ...
                                 'posixtime')),29);
if strcmp(tstr1,tstr2)
    tstr = tstr1;
else
    tstr = [tstr1,' -- ',tstr2];
end
% title(h1,tstr)
drawnow


end