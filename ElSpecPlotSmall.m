function fignum = ElSpecPlotSmall(ElSpecOut,varargin)
%
% Plot ElSpec Output, this function plots only the energy spectrum,
% FAC, and auroral power
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
%    faclim    FAC limits (mA/m^2), default [0 10]
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

defaultElim = [1 100];%[min(ElSpecOut.Ec) max(ElSpecOut.Ec)]./1000;
checkElim = @(x) (isnumeric(x) & (length(x)==2));

defaultFAClim = [0 10];
checkFAClim = @(x) (isnumeric(x) & (length(x)==2));

defaultPlim = [0 50];
checkPlim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecOut.te(1) + ElSpecOut.ts(1))/2,'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecOut.te(end) + ElSpecOut.ts(end))/2,'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultFluxtype = 'energy';
validFluxtype = {'number','energy','both'};
checkFluxtype = @(x) any(validatestring(x,validFluxtype));

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
addParameter(p,'faclim',defaultFAClim,checkFAClim)
addParameter(p,'plim',defaultPlim,checkPlim)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'fluxtype',defaultFluxtype,checkFluxtype)
addParameter(p,'emin',defaultEmin,checkEmin)
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
% heighs
hh = ElSpecOut.h;
% energies
EE = ElSpecOut.Ec/1000;

% update FAC and power using the Emin input
if isfield(ElSpecOut,'Emin')
    emin0 = ElSpecOut.Emin;
elseif isfield(ElSpecOut,'emin')
    emin0 = ElSpecOut.emin;
else
    emin0 = 0;
end
if p.Results.emin~=emin0
    Eind_fac = ElSpecOut.Ec >= p.Results.emin;
    for tind=1:length(ElSpecOut.ts)
        ElSpecOut.FAC(tind) = sum(ElSpecOut.Ie(Eind_fac,tind).* ...
                                  ElSpecOut.dE(Eind_fac)')*1.60217662e-19;
        if isfield(ElSpecOut,'IeCov')
            ElSpecOut.FACstd(tind) = sqrt(ElSpecOut.dE(Eind_fac) * ...
                                          ElSpecOut.IeCov(Eind_fac,Eind_fac,tind) * ElSpecOut.dE(Eind_fac)')*1.60217662e-19;
        else
            ElSpecOut.FACstd(tind) = 0;
        end

        % Power carried by the precipitating electrons
        ElSpecOut.Pe(tind) = sum(ElSpecOut.Ie(Eind_fac,tind).*ElSpecOut.dE(Eind_fac)'.*ElSpecOut.Ec(Eind_fac)')*1.60217662e-19;
        EdE = ElSpecOut.dE(Eind_fac)'.*ElSpecOut.Ec(Eind_fac)';
        if isfield(ElSpecOut,'IeCov')
            ElSpecOut.PeStd(tind) = sqrt( EdE' * ...
                                          ElSpecOut.IeCov(Eind_fac,Eind_fac,tind)* EdE ) * 1.60217662e-19;
        else
            ElSpecOut.PeStd(tind) = 0;
        end

    end
end

% number vs energy flux
switch lower(p.Results.fluxtype)

  case 'number'
    Ie = ElSpecOut.Ie;
    Ielabel={'[(eVm^2s)^{-1}]'};
    ielim = p.Results.ielim;
  case 'energy'
    Ie = ElSpecOut.Ie;
    for iiE = 1:length(ElSpecOut.Ec)
        Ie(iiE,:) = Ie(iiE,:).*ElSpecOut.Ec(iiE);
    end
    Ielabel={'[eV/(eVm^2s)]'};
    ielim = p.Results.ieelim;
  case 'both'
    Ie = ElSpecOut.Ie;
    Ielabel={'[(eVm^2s)^{-1}]'};
    IEe = ElSpecOut.Ie;
    for iiE = 1:length(ElSpecOut.Ec)
        IEe(iiE,:) = IEe(iiE,:).*ElSpecOut.Ec(iiE);
    end
    IEelabel={'[eV/(eVm^2s)]'};
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
colormap(jet);

h1=subplot(3,1,1);
pcolor(ts,EE,log10(Ie)),shading flat
set(h1,'yscale','log')
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
caxis(ielim)
ylabel('Energy [keV]')
cbh1=colorbar;
ylabel(cbh1,Ielabel)

h2=subplot(3,1,2);
plot(tt,ElSpecOut.FAC*NaN)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.faclim)
hold on
for tind=1:length(ElSpecOut.ts)
    plot([tt(tind) tt(tind)] , [-1 1]*ElSpecOut.FACstd(tind)*1e6+ElSpecOut.FAC(tind)*1e6,'r-','LineWidth',1)
end
plot(tt,ElSpecOut.FAC*1e6,'k-','LineWidth',1)
hold off
ylabel('FAC [\muAm^{-2}]')
grid on


h3=subplot(3,1,3);
plot(tt,ElSpecOut.Pe*NaN)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.plim)
hold on
for tind=1:length(ElSpecOut.ts)
    plot([tt(tind) tt(tind)] , [-1 1]*ElSpecOut.PeStd(tind)*1e3+ElSpecOut.Pe(tind)*1e3,'r-','LineWidth',1)
end
plot(tt,ElSpecOut.Pe*1000,'k-','LineWidth',1)
hold off
ylabel('Power [mWm^{-2}]')
xlabel('Time [UTC]');
grid on

drawnow % this is important to update everything before setting the sizes...
pos1=get(h1,'Position');
pos2=get(h2,'Position');
pos3=get(h3,'Position');
set(h1, 'Position', [pos1(1:2)-[0,.06] [.85 1.4].*pos3(3:4)]);
set(h2, 'Position', [pos2(1:2)-[0,.04] [.85 1.4].*pos3(3:4)]);
set(h3, 'Position', [pos3(1:2)-[0,.02] [.85 1.4].*pos3(3:4)]);
linkaxes([h1 h2 h3],'x')

eticks=[.5 1 2 4 8 16 32 64 128 256 512 1024];
set(h1,'YTick',eticks(eticks<p.Results.elim(2)))
set(h1,'TickDir','both')
set(h2,'TickDir','both')
facstep = 1;
if range(p.Results.faclim) > 4
    facstep = 2;
end
if range(p.Results.faclim) > 10
    facstep = 4;
end
if range(p.Results.faclim) > 15
    facstep = 5;
end
if range(p.Results.faclim) > 30
    facstep = 10;
end
facticks = ceil(p.Results.faclim(1)):facstep:floor(p.Results.faclim(end)-.1);
set(h2,'YTick',facticks)
pstep = 1;
if range(p.Results.plim) > 4
    pstep = 2;
end
if range(p.Results.plim) > 10
    pstep = 4;
end
if range(p.Results.plim) > 15
    pstep = 5;
end
if range(p.Results.plim) > 30
    pstep = 10;
end
if range(p.Results.plim) > 60
    pstep = 20;
end
if range(p.Results.plim) > 150
    pstep = 50;
end
pticks = ceil(p.Results.plim(1)):pstep:floor(p.Results.plim(end)-.1);
set(h3,'YTick',pticks)
set(h3,'TickDir','both')

set(h2,'XTick',get(h1,'XTick'))
set(h3,'XTick',get(h1,'XTick'))
datetick(h2,'x',13,'keeplimits','keepticks')
datetick(h3,'x',13,'keeplimits','keepticks')
set(h1,'XTickLabel','')
set(h2,'XTickLabel','')
set(h1,'LineWidth',1)
set(h2,'LineWidth',1)
set(h3,'LineWidth',1)

cbsize1=get(cbh1,'Position');
set(cbh1,'Position',cbsize1.*[1,1.02,1,.8]);

set(h1,'fontsize',12)
set(cbh1,'fontsize',12)
set(h2,'fontsize',12)
set(h3,'fontsize',12)

set(cbh1,'YTick',5:13)
set(cbh1,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})


set(h1,'layer','top')
set(cbh1,'LineWidth',1,'TickDir','both')
set(h2,'layer','top')
set(h3,'layer','top')

tstr1 = datestr(datenum(datetime((ElSpecOut.te(1)),'ConvertFrom', ...
                                 'posixtime')),29);
tstr2 = datestr(datenum(datetime((ElSpecOut.te(end)),'ConvertFrom', ...
                                 'posixtime')),29);
if strcmp(tstr1,tstr2)
    tstr = tstr1;
else
    tstr = [tstr1,' -- ',tstr2];
end
title(h1,tstr)
drawnow


end