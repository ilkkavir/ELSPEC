function fignum = ElSpecPlot(ElSpecOut,varargin)
%
% Plot ElSpec Output.
%
% fignum = ElSpecPlot( ElSpecOut , <fignum> , ... )
%
% INPUT:
%  ElSpecOut   an output list from ElSpec
%  fignum      optional number of the figure. A new figure is
%              created by default.
%  ...         optional parameters as name-value pairs:
%    nelim     electron density color scale limits, (log(m^-3)), default [10 12]
%    ielim     differential number flux color scale limits, (log((eVm^2s)^-1)),
%              default [6 10]
%    ieelim    differential energy flux color scale limits
%              (log(eV/(eVm^2s))), default [10 14]
%    elim      energy axis limits for the differential flux plot
%              (keV), default [1 100]
%    faclim    FAC limits (mA/m^2), default [0 10]
%    plim      net energy flux limits (W/^2), default [0 50]
%    chisqrlim chi-square plot limits, default [0 10]
%    btime     plot start time as vector [yyyy mm dd HH MM SS],
%              default: first data point
%    etime     plot end time as vector [yyyy mm dd HH MM SS],
%              default: last data point
%    fontsize  font size in axis labels. Default 8
%    fluxtype  flux plot type, 'number', 'energy', or 'both'. If
%              'both' is selected, the net energy flux plot is
%              omitted. Default: 'energy'
%    neplot    type of the electron density plot, 'log' or
%              'linear'. default 'log'
%    emin      minimum energy to include in the FAC and net energy
%              flux integrations (eV). Default 1000. NOTICE: error
%              estimates cannot be calculated if 'emin' is changed
%              and ElSpecOut was calculated with
%              'saveiecov','false' in ElSpec.
%    cutgaps   logical, should white space be plotted on obvious
%              data gaps. Default true
%    visible   figure visibility. Use 'on' to plot on a visible device,
%              'off' to avoid opening a figure window. Default 'on'
%    plimE0    The peak energy is plotted as a line in the spectrum panel,
%              when the integral energy flux is larger then plimE0 mW/m^2
%
%
% IV 2016, 2017, 2018
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


p = inputParser;

defaultFignum = -1;

defaultHlim = [ElSpecOut.h(1) ElSpecOut.h(end)];
checkHlim = @(x) (isnumeric(x) & (length(x)==2));

defaultNelim = [10 12];
checkNelim = @(x) (isnumeric(x) & (length(x)==2));

defaultIelim = [6 10];
checkIelim = @(x) (isnumeric(x) & (length(x)==2));

defaultIEelim = [10 14];
checkIEelim = @(x) (isnumeric(x) & (length(x)==2));

defaultElim = [1 100];%[min(ElSpecOut.Ec) max(ElSpecOut.Ec)]./1000;
checkElim = @(x) (isnumeric(x) & (length(x)==2));

%defaultFAClim = [0 10];
%defaultFAClim = [0 ceil(max(ElSpecOut.FAC)*2e5)*5];
defaultFAClim = [0 ceil( min( max(ElSpecOut.FAC) , quantile(ElSpecOut.FAC,.99)*2) * 2e5)*5];
checkFAClim = @(x) (isnumeric(x) & (length(x)==2));

%defaultPlim = [0 50];
%defaultPlim = [0 ceil(max(ElSpecOut.Pe)*50)*20];
defaultPlim = [0 ceil( min( max(ElSpecOut.Pe) , quantile(ElSpecOut.Pe,.99)*2) * 50)*20];
checkPlim = @(x) (isnumeric(x) & (length(x)==2));

defaultChisqrlim = [0 10];
checkChisqrlim = @(x) (isnumeric(x) & (length(x)==2));

defaultFontsize = 8;
checkFontsize = @(x) (isnumeric(x) & (length(x)==1));

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

defaultVisible = 'on';
checkVisible = @(x) (ischar(x)); % this is enough,the option will be checked by figure()

defaultPlimE0 = 2;
checkPlimE0 = @(x) (isnumeric(x));

if isfield(ElSpecOut,'Emin')
    defaultEmin = ElSpecOut.Emin;
elseif isfield(ElSpecOut,'emin')
    defaultEmin = ElSpecOut.emin;
else
    defaultEmin = 1e3;
end
checkEmin = @(x) (isnumeric(x) & length(x)==1 );

addRequired(p,'ElSpecOut',@isstruct);
addOptional(p,'fignum',defaultFignum,@isnumeric);
addParameter(p,'hlim',defaultHlim,checkHlim)
addParameter(p,'nelim',defaultNelim,checkNelim)
addParameter(p,'ielim',defaultIelim,checkIelim)
addParameter(p,'ieelim',defaultIEelim,checkIEelim)
addParameter(p,'elim',defaultElim,checkElim)
addParameter(p,'faclim',defaultFAClim,checkFAClim)
addParameter(p,'plim',defaultPlim,checkPlim)
addParameter(p,'chisqrlim',defaultChisqrlim,checkChisqrlim)
addParameter(p,'fontsize',defaultFontsize,checkFontsize)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'fluxtype',defaultFluxtype,checkFluxtype)
addParameter(p,'neplot',defaultNeplot,checkNeplot)
addParameter(p,'emin',defaultEmin,checkEmin)
addParameter(p,'cutgaps',defaultCutgaps,checkCutgaps)
addParameter(p,'visible',defaultVisible,checkVisible)
addParameter(p,'plime0',defaultPlimE0,checkPlimE0)
parse(p,ElSpecOut,varargin{:})


% put white space on data gaps
if p.Results.cutgaps
    mediandt = median(diff(ElSpecOut.ts));
    gaps = diff(ElSpecOut.ts) > (5*mediandt);
    rminds = [find(gaps) ; (find(gaps)+1)];
    ElSpecOut.pp(:,rminds) = NaN;
    ElSpecOut.ne(:,rminds) = NaN;
    ElSpecOut.Ie(:,rminds) = NaN;
    ElSpecOut.FAC(rminds) = NaN;
    ElSpecOut.FACstd(rminds) = NaN;
    ElSpecOut.Pe(rminds) = NaN;
    ElSpecOut.PeStd(rminds) = NaN;
    ElSpecOut.chisqr(rminds) = NaN;
    ElSpecOut.E0(rminds) = NaN;
end


% centres of time-bins as matlab datenum
tt = datenum(datetime((ElSpecOut.te + ElSpecOut.ts)/2,'ConvertFrom','posixtime'));
% start points of time-bins for the pcolor plots
ts = datenum(datetime(ElSpecOut.ts,'ConvertFrom','posixtime'));
% end points of time-bins
te = datenum(datetime(ElSpecOut.te,'ConvertFrom','posixtime'));
% heighs
hh = ElSpecOut.h;
% energies
EE = ElSpecOut.Ec/1000;
% power profiles with negatives removed
pp = ElSpecOut.pp;
pp(pp<1) = 1;

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

% for testing that what fluxes cannot be distinguished from zero
%naninds = ElSpecOut.IeStd > ElSpecOut.Ie;
%Ie(naninds) = NaN;



fignum = p.Results.fignum;

if fignum>0
    fh = figure(fignum);
    figpos = get(fh,'Position');
else
    fh = figure('visible',p.Results.visible);
    fignum = fh.Number;
    figpos = [0 0 21/29.7*690 690];
end
figpos(3:4) = [21/29.7*690 690];
set(fh,'Position',figpos,'PaperPositionMode','Auto');
colormap(jet);

h1=subplot(6,1,1);
if strcmp(p.Results.neplot,'log')
    pcolor(ts,hh,log10(pp)),shading flat
    caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
    pcolor(ts,hh,pp),shading flat
    caxis(10.^(p.Results.nelim))
else
error(['Unknown neplot ',neplot])
end
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%datetick('x',13,'keeplimits')
ylabel('Height [km]')
cbh1=colorbar;
ylabel(cbh1,{'N_e [m^{-3}] (meas.)'})

h2=subplot(6,1,2);
if strcmp(p.Results.neplot,'log')
    ppOut = real(ElSpecOut.ne); % Ne should be real, but we have seen imaginary parts with unrealistic inputs...
    ppOut(ppOut<1) = 1;
    pcolor(ts,hh,log10(ppOut)),shading flat
    caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
    pcolor(ts,hh,real(ElSpecOut.ne)),shading flat
    caxis(10.^(p.Results.nelim))
else
    error(['Unknown neplot ',neplot])
end
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.hlim)
cbh2=colorbar;
ylabel('Height [km]')
ylabel(cbh2,{'N_e [m^{-3}] (model)'})

h3=subplot(6,1,3);
pcolor(ts,EE,log10(Ie)),shading flat
set(h3,'yscale','log')
if isfinite(p.Results.plime0)
    hold on
    E0 = ElSpecOut.E0/1e3;
    E0(ElSpecOut.Pe*1000<p.Results.plime0) = NaN;
    %    E0(ElSpecOut.Pe*1e3 < .01) = NaN;
    %    E0 = [E0(:),E0(:)]';
    E0 = E0(:);
    %    tttmp = [ts(:),ts(:)]';
    %    tttmp = tttmp(:);
    %    ttE0 = [tttmp(2:end);te(end)];
    %    plot(ttE0,E0,'m-','linewidth',.5)
    plot(tt,E0,'m.')
    hold off
end
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(ielim)
ylabel('Energy [keV]')
cbh3=colorbar;
ylabel(cbh3,Ielabel)
%ylim([.5 20])

h4=subplot(6,1,4);
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


h5=subplot(6,1,5);
if strcmp(p.Results.fluxtype,'both')
    pcolor(ts,EE,log10(IEe)),shading flat
    set(h5,'yscale','log')
    ylim(p.Results.elim)
    xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
    caxis(p.Results.ieelim)
    ylabel('Energy [keV]')
    cbh5=colorbar;
    ylabel(cbh5,IEelabel)
else
    plot(tt,ElSpecOut.Pe*NaN)
    xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
    ylim(p.Results.plim)
    hold on
    for tind=1:length(ElSpecOut.ts)
        plot([tt(tind) tt(tind)] , [-1 1]*ElSpecOut.PeStd(tind)*1e3+ElSpecOut.Pe(tind)*1e3,'r-','LineWidth',1)
    end
    plot(tt,ElSpecOut.Pe*1000,'k-','LineWidth',1)
    hold off
    %    ylabel('Power [mWm^{-2}]')
    ylabel({'Energy flux';'[mWm^{-2}]'})
end
grid on

h6=subplot(6,1,6);
plot(tt,ElSpecOut.chisqr,'k-','LineWidth',1)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.chisqrlim);
ylabel('\chi^2')
xlabel('Time [UTC]');
grid on


drawnow % this is important to update everything before setting the sizes...


xlims = datenum([datetime(p.Results.btime) datetime(p.Results.etime)]);
xstep = 1; % time step in seconds
nx = diff(xlims)/(xstep/3600/24);
tickfmt = 'HH:MM:SS';
while nx > 10
    if xstep == 1; xstepn = 2; end
    if xstep == 2; xstepn = 5; end
    if xstep == 5; xstepn = 10; end
    if xstep == 10; xstepn = 15; end
    if xstep == 15; xstepn = 30; end
    if xstep == 30; xstepn = 60; end
    if xstep == 60; xstepn = 2*60; end
    if xstep == 2*60; xstepn = 5*60; end
    if xstep == 5*60; xstepn = 10*60; end
    if xstep == 10*60; xstepn = 15*60; end
    if xstep == 15*60; xstepn = 30*60; end
    if xstep == 30*60; xstepn = 60*60; end
    if xstep == 60*60; xstepn = 2*60*60; end
    if xstep == 2*60*60; xstepn = 3*60*60; end
    if xstep == 3*60*60; xstepn = 4*60*60; end
    if xstep == 4*60*60; xstepn = 5*60*60; end
    if xstep == 5*60*60; xstepn = 6*60*60; end
    if xstep == 6*60*60; xstepn = 12*60*60; end
    if xstep == 12*60*60; xstepn = 24*60*60; end
    xstep = xstepn;
    if xstep >= 60
        tickfmt = "HH:MM";
    end
    if xstep >= 3600
        tickfmt = "HH";
    end
    nx = diff(xlims)/(xstep/3600/24);
end
xstep = xstep/3600/24;
xticks = round(xlims(1)/xstep)*xstep : xstep : round(xlims(2)/xstep)*xstep;
xticks = xticks(xticks>=xlims(1));
xticks = xticks(xticks<=xlims(2));
set(h1,'xtick',xticks);
set(h2,'xtick',xticks);
set(h3,'xtick',xticks);
set(h4,'xtick',xticks);
set(h5,'xtick',xticks);
set(h6,'xtick',xticks);

pos1=get(h1,'Position');
pos2=get(h2,'Position');
pos3=get(h3,'Position');
pos4=get(h4,'Position');
pos5=get(h5,'Position');
pos6=get(h6,'Position');
panelshift = -[0,.05];
panelscale = [.85 1.65];
set(h1, 'Position', [pos1(1:2)+panelshift panelscale.*pos4(3:4)]);
set(h2, 'Position', [pos2(1:2)+panelshift panelscale.*pos4(3:4)]);
set(h3, 'Position', [pos3(1:2)+panelshift panelscale.*pos4(3:4)]);
set(h4, 'Position', [pos4(1:2)+panelshift panelscale.*pos4(3:4)]);
set(h5, 'Position', [pos5(1:2)+panelshift panelscale.*pos4(3:4)]);
set(h6, 'Position', [pos6(1:2)+panelshift panelscale.*pos4(3:4)]);
linkaxes([h1 h2 h3 h4 h5 h6],'x')
linkaxes([h1 h2],'y')
if strcmp(p.Results.fluxtype,'both')
    linkaxes([h3 h5],'y')
end

hstep = 1;
if range(p.Results.hlim)>2
    hstep = 2;
end
if range(p.Results.hlim)>10
    hstep = 5;
end
if range(p.Results.hlim)>20
    hstep = 10;
end
if range(p.Results.hlim)>50
    hstep = 20;
end
hticks = ceil(p.Results.hlim(1)/10)*10:hstep:floor(p.Results.hlim(2)/10-.01)*10;
set(h1,'TickDir','both')
set(h1,'YTick',hticks)
set(h2,'TickDir','both')
set(h2,'YTick',hticks)
eticks=[.5 1 2 4 8 16 32 64 128 256 512 1024];
set(h3,'YTick',eticks(eticks<p.Results.elim(2)))
set(h3,'TickDir','both')
set(h4,'TickDir','both')
facstep = .1;
if range(p.Results.faclim) > .5
    facstep = .2;
end
if range(p.Results.faclim) > 1
    facstep = .5;
end
if range(p.Results.faclim) > 2.5
    facstep = 1;
end
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
facticks = ceil(p.Results.faclim(1)*10)/10:facstep:floor(p.Results.faclim(end)*10-.1)/10;
if facstep >=2
    facticks = ceil(p.Results.faclim(1)):facstep:floor(p.Results.faclim(end)-.1);
end
set(h4,'YTick',facticks)
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
if strcmp(p.Results.fluxtype,'both')
    set(h5,'YTick',eticks(eticks<p.Results.elim(2)))
else
    set(h5,'YTick',pticks)
end
set(h5,'TickDir','both')
chistep=1;
if range(p.Results.chisqrlim) > 5
    chistep=2;
end
if range(p.Results.chisqrlim) > 10
    chistep=4;
end
chiticks = ceil(p.Results.chisqrlim(1)):chistep:floor(p.Results.chisqrlim(end)-.1);
set(h6,'TickDir','both')
set(h6,'Ytick',chiticks)
%set(h2,'XTick',get(h1,'XTick'))
%set(h3,'XTick',get(h1,'XTick'))
%set(h4,'XTick',get(h1,'XTick'))
%set(h5,'XTick',get(h1,'XTick'))
%set(h6,'XTick',get(h1,'XTick'))
% datetick(h1,'x',13,'keeplimits')
% datetick(h2,'x',13,'keeplimits')
% datetick(h3,'x',13,'keeplimits')
% datetick(h4,'x',13,'keeplimits')
% datetick(h5,'x',13,'keeplimits')
% datetick(h6,'x',13,'keeplimits')
datetick(h1,'x',tickfmt,'keeplimits','keepticks')
datetick(h2,'x',tickfmt,'keeplimits','keepticks')
datetick(h3,'x',tickfmt,'keeplimits','keepticks')
datetick(h4,'x',tickfmt,'keeplimits','keepticks')
datetick(h5,'x',tickfmt,'keeplimits','keepticks')
datetick(h6,'x',tickfmt,'keeplimits','keepticks')
set(h1,'XTickLabel','')
set(h2,'XTickLabel','')
set(h3,'XTickLabel','')
set(h4,'XTickLabel','')
set(h5,'XTickLabel','')
set(h1,'LineWidth',1)
set(h2,'LineWidth',1)
set(h3,'LineWidth',1)
set(h4,'LineWidth',1)
set(h5,'LineWidth',1)
set(h6,'LineWidth',1)

cbsize1=get(cbh1,'Position');
set(cbh1,'Position',cbsize1.*[1,1.02,1,.8]);

cbsize2=get(cbh2,'Position');
set(cbh2,'Position',cbsize2.*[1,1.02,1,.8]);

cbsize3=get(cbh3,'Position');
set(cbh3,'Position',cbsize3.*[1,1.02,1,.8]);


set(h1,'fontsize',p.Results.fontsize)
set(cbh1,'fontsize',p.Results.fontsize)
set(h2,'fontsize',p.Results.fontsize)
set(cbh2,'fontsize',p.Results.fontsize)
set(h3,'fontsize',p.Results.fontsize)
set(cbh3,'fontsize',p.Results.fontsize)
set(h4,'fontsize',p.Results.fontsize)
set(h5,'fontsize',p.Results.fontsize)
set(h6,'fontsize',p.Results.fontsize)

if strcmp(p.Results.neplot,'log')
    set(cbh1,'YTick',5:13)
    set(cbh1,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})
    set(cbh2,'YTick',5:13)
    set(cbh2,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'}) ...

end
set(cbh3,'YTick',5:13)
set(cbh3,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})


set(h1,'layer','top')
set(cbh1,'LineWidth',1,'TickDir','both')
set(h2,'layer','top')
set(cbh2,'LineWidth',1,'TickDir','both')
set(h3,'layer','top')
set(cbh3,'LineWidth',1,'TickDir','both')
set(h4,'layer','top')
set(h5,'layer','top')
set(h6,'layer','top')


% time string corresponding to start and end of the measurement. One percent tolerance to
% avoid double time stamps in cases where just one time step is from the previous/next day
tdiff = ElSpecOut.te(end)-ElSpecOut.te(1);

tstr1 = datestr(datenum(datetime((ElSpecOut.te(1)+.01*tdiff),'ConvertFrom', ...
                                 'posixtime')),29);
tstr2 = datestr(datenum(datetime((ElSpecOut.te(end)-.01*tdiff),'ConvertFrom', ...
                                 'posixtime')),29);


if strcmp(tstr1,tstr2)
    tstr = tstr1;
else
    tstr = [tstr1,' -- ',tstr2];
end
%title(h1,tstr)
title(h1,[upper(ElSpecOut.radar) ' ' tstr])
drawnow


end


function rr = range(x)
    rr = max(x) - min(x);
end
