function fignum = ElSpecGlowPlot(ElSpecGlowOut,varargin)
%
% Plot ElSpec and glow emission profile results
%
% fignum = ElSpecGlowPlot( ElSpecGlowOut , <fignum> , ... )
%
% INPUT:
%  ElSpecGlowOut   an output list from GLOWemissions
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
%    plim      net energy flux limits (W/^2), default [0 50]
%    e5577lim  color scale limits for the 5577 Å emission, default [0 3e4]
%    btime     plot start time as vector [yyyy mm dd HH MM SS],
%              default: first data point
%    etime     plot end time as vector [yyyy mm dd HH MM SS],
%              default: last data point
%    fluxtype  flux plot type, 'number' or 'energy'. Default: 'energy'
%    neplot    type of the electron density plot, 'log' or
%              'linear'. default 'log'
%    emin      minimum energy to include in the  net energy
%              flux integrations (eV). Default 1000. NOTICE: error
%              estimates cannot be calculated if 'emin' is changed
%              and ElSpecGlowOut was calculated with
%              'saveiecov','false' in ElSpec.
%    cutgaps   logical, should white space be plotted on obvious
%              data gaps. Default true
%
%
% IV 2016, 2017, 2018
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


p = inputParser;

defaultFignum = -1;

defaultHlim = [ElSpecGlowOut.h(1) ElSpecGlowOut.h(end)];
checkHlim = @(x) (isnumeric(x) & (length(x)==2));

defaultNelim = [10 12];
checkNelim = @(x) (isnumeric(x) & (length(x)==2));

defaultIelim = [6 10];
checkIelim = @(x) (isnumeric(x) & (length(x)==2));

defaultIEelim = [10 14];
checkIEelim = @(x) (isnumeric(x) & (length(x)==2));

defaultElim = [1 100];%[min(ElSpecGlowOut.Ec) max(ElSpecGlowOut.Ec)]./1000;
checkElim = @(x) (isnumeric(x) & (length(x)==2));

defaultPlim = [0 50];
checkPlim = @(x) (isnumeric(x) & (length(x)==2));

defaultChisqrlim = [0 10];
checkChisqrlim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecGlowOut.te(1) + ElSpecGlowOut.ts(1))/2,'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecGlowOut.te(end) + ElSpecGlowOut.ts(end))/2,'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultFluxtype = 'energy';
validFluxtype = {'number','energy'};
checkFluxtype = @(x) any(validatestring(x,validFluxtype));

defaultNeplot = 'log';
validNeplot = {'log','linear'};
checkNeplot = @(x) any(validatestring(x,validNeplot));

defaultCutgaps = 1;
checkCutgaps = @(x) (islogical(x)|isnumeric(x));

defaultE5577lim = [0 3e4];
checkE5577lim = @(x) (isnumeric(x) & length(x)==2);

if isfield(ElSpecGlowOut,'Emin')
    defaultEmin = ElSpecGlowOut.Emin;
else
    defaultEmin = 1e3;
end
checkEmin = @(x) (isnumeric(x) & length(x)==1 );

addRequired(p,'ElSpecGlowOut',@isstruct);
addOptional(p,'fignum',defaultFignum,@isnumeric);
addParameter(p,'hlim',defaultHlim,checkHlim)
addParameter(p,'nelim',defaultNelim,checkNelim)
addParameter(p,'ielim',defaultIelim,checkIelim)
addParameter(p,'ieelim',defaultIEelim,checkIEelim)
addParameter(p,'elim',defaultElim,checkElim)
addParameter(p,'plim',defaultPlim,checkPlim)
addParameter(p,'chisqrlim',defaultChisqrlim,checkChisqrlim)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'fluxtype',defaultFluxtype,checkFluxtype)
addParameter(p,'neplot',defaultNeplot,checkNeplot)
addParameter(p,'emin',defaultEmin,checkEmin)
addParameter(p,'cutgaps',defaultCutgaps,checkCutgaps)
addParameter(p,'e5577lim',defaultE5577lim,checkE5577lim)
parse(p,ElSpecGlowOut,varargin{:})


% put white space on data gaps
if p.Results.cutgaps
    mediandt = median(diff(ElSpecGlowOut.ts));
    gaps = diff(ElSpecGlowOut.ts) > (10*mediandt);
    rminds = [find(gaps) ; (find(gaps)+1)];
    ElSpecGlowOut.pp(:,rminds) = NaN;
    ElSpecGlowOut.ne(:,rminds) = NaN;
    ElSpecGlowOut.Ie(:,rminds) = NaN;
    ElSpecGlowOut.Pe(rminds) = NaN;
    ElSpecGlowOut.PeStd(rminds) = NaN;
    ElSpecGlowOut.chisqr(rminds) = NaN;
    ElSpecGlowOut.GLOWne(:,rminds) = NaN;
    ElSpecGlowOut.GLOWe5577(:,rminds) = NaN;
    ElSpecGlowOut.GLOWe5577max(rminds) = NaN;
    ElSpecGlowOut.GLOWh5577max(rminds) = NaN;
    ElSpecGlowOut.GLOWh5577mean(rminds) = NaN;
end


% centres of time-bins as matlab datenum
tt = datenum(datetime((ElSpecGlowOut.te + ElSpecGlowOut.ts)/2,'ConvertFrom','posixtime'));
% start points of time-bins for the pcolor plots
ts = datenum(datetime(ElSpecGlowOut.ts,'ConvertFrom','posixtime'));
% heighs
hh = ElSpecGlowOut.h;
% energies
EE = ElSpecGlowOut.Ec/1000;
% power profiles with negatives removed
pp = ElSpecGlowOut.pp;
pp(pp<1) = 1;

% update  power using the Emin input
if isfield(ElSpecGlowOut,'Emin')
    emin0 = ElSpecGlowOut.Emin;
elseif isfield(ElSpecGlowOut,'emin')
    emin0 = ElSpecGlowOut.emin;
else
    emin0 = 0;
end
if p.Results.emin~=emin0
    Eind_fac = ElSpecGlowOut.Ec >= p.Results.emin;
    for tind=1:length(ElSpecGlowOut.ts)
        % Power carried by the precipitating electrons
        ElSpecGlowOut.Pe(tind) = sum(ElSpecGlowOut.Ie(Eind_fac,tind).*ElSpecGlowOut.dE(Eind_fac)'.*ElSpecGlowOut.Ec(Eind_fac)')*1.60217662e-19;
        EdE = ElSpecGlowOut.dE(Eind_fac)'.*ElSpecGlowOut.Ec(Eind_fac)';
        if isfield(ElSpecGlowOut,'IeCov')
            ElSpecGlowOut.PeStd(tind) = sqrt( EdE' * ...
                                          ElSpecGlowOut.IeCov(Eind_fac,Eind_fac,tind)* EdE ) * 1.60217662e-19;
        else
            ElSpecGlowOut.PeStd(tind) = 0;
        end

    end
end

% number vs energy flux
switch lower(p.Results.fluxtype)

  case 'number'
    Ie = ElSpecGlowOut.Ie;
    Ielabel={'[(eVm^2s)^{-1}]'};
    ielim = p.Results.ielim;
  case 'energy'
    Ie = ElSpecGlowOut.Ie;
    for iiE = 1:length(ElSpecGlowOut.Ec)
        Ie(iiE,:) = Ie(iiE,:).*ElSpecGlowOut.Ec(iiE);
    end
    Ielabel={'[eV/(eVm^2s)]'};
    ielim = p.Results.ieelim;
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
figpos(3:4) = [21/29.7*690 690];
set(gcf,'Position',figpos,'PaperPositionMode','Auto','color','white');
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
datetick('x',13,'keeplimits')
yh1=ylabel('Height [km]');
cbh1=colorbar;
ycbh1=ylabel(cbh1,{'N_e [m^{-3}] (meas.)'});

h2=subplot(6,1,2);
if strcmp(p.Results.neplot,'log')
    ppOut = ElSpecGlowOut.ne;
    ppOut(ppOut<1) = 1;
    pcolor(ts,hh,log10(ppOut)),shading flat
    caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
    pcolor(ts,hh,ElSpecGlowOut.ne),shading flat
    caxis(10.^(p.Results.nelim))
else
    error(['Unknown neplot ',neplot])
end
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.hlim)
cbh2=colorbar;
yh2=ylabel('Height [km]');
ycbh2=ylabel(cbh2,{'N_e [m^{-3}] (ELSPEC)'});


% GLOW Ne
h3=subplot(6,1,3);
if strcmp(p.Results.neplot,'log')
    pcolor(ts,ElSpecGlowOut.GLOWh,log10(ElSpecGlowOut.GLOWne)),shading flat
    caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
    pcolor(ts,ElSpecGlowOut.GLOWh,ElSpecGlowOut.GLOWne),shading flat
    caxis(10.^(p.Results.nelim))
else
error(['Unknown neplot ',neplot])
end
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
yh3=ylabel('Height [km]');
cbh3=colorbar;
ycbh3=ylabel(cbh3,{'N_e [m^{-3}] (GLOW)'});


% energy spectrum
h4=subplot(6,1,4);
pcolor(ts,EE,log10(Ie)),shading flat
set(h4,'yscale','log')
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(ielim)
yh4=ylabel('Energy [keV]');
cbh4=colorbar;
ycbh4=ylabel(cbh4,Ielabel);
%ylim([.5 20])

% emission rate from GLOW
h5=subplot(6,1,5);
pcolor(ts,ElSpecGlowOut.GLOWh,ElSpecGlowOut.GLOWe5577),shading flat
caxis(p.Results.e5577lim)
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
datetick('x',13,'keeplimits')
yh5=ylabel('Height [km]');
cbh5=colorbar;
ycbh5=ylabel(cbh5,{'557.7 nm [cm^{-3}s^{-1}]'});
hold on
l1=plot(datenum(ElSpecGlowOut.GLOWtime),ElSpecGlowOut.GLOWh5577max,'r-','linewidth',1.5);
%l2=plot(datenum(ElSpecGlowOut.GLOWtime),ElSpecGlowOut.GLOWh5577mean,'w-','linewidth',1.5);
%legend([l1,l2],{'Peak','Mean'},'location','southwest','textcolor','w')
legend([l1],{'Peak'},'location','southwest','textcolor','w')
legend('boxoff')


% auroral power from ELSPEC
h6=subplot(6,1,6);
plot(tt,ElSpecGlowOut.Pe*NaN)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.plim)
hold on
for tind=1:length(ElSpecGlowOut.ts)
    plot([tt(tind) tt(tind)] , [-1 1]*ElSpecGlowOut.PeStd(tind)*1e3+ElSpecGlowOut.Pe(tind)*1e3,'r-','LineWidth',1)
end
plot(tt,ElSpecGlowOut.Pe*1000,'k-','LineWidth',1)
hold off
xlabel('Time [UTC]');
yh6=ylabel('Power [mWm^{-2}]');
grid on



drawnow % this is important to update everything before setting the sizes...
pos1=get(h1,'Position');
pos2=get(h2,'Position');
pos3=get(h3,'Position');
pos4=get(h4,'Position');
pos5=get(h5,'Position');
pos6=get(h6,'Position');
set(h1, 'Position', [pos1(1:2)-[0,.05] 0.66518  0.12874]);
set(h2, 'Position', [pos2(1:2)-[0,.05] 0.66518  0.12874]);
set(h3, 'Position', [pos3(1:2)-[0,.05] 0.66518  0.12874]);
set(h4, 'Position', [pos4(1:2)-[0,.05] 0.66518  0.12874]);
set(h5, 'Position', [pos5(1:2)-[0,.05] 0.66518  0.12874]);
set(h6, 'Position', [pos6(1:2)-[0,.05] 0.66518  0.12874]);

drawnow

linkaxes([h1 h2 h3 h4 h5 h6],'x')
linkaxes([h1 h2,h3,h5],'y')

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
set(h3,'YTick',hticks)
set(h3,'TickDir','both')
eticks=[.5 1 2 4 8 16 32 64 128 256 512 1024];
set(h4,'YTick',eticks(eticks<p.Results.elim(2)))
set(h4,'TickDir','both')
set(h5,'YTick',hticks)
set(h5,'TickDir','both')
set(h6,'TickDir','both')

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
set(h6,'YTick',pticks)

drawnow

set(h2,'XTick',get(h1,'XTick'))
set(h3,'XTick',get(h1,'XTick'))
set(h4,'XTick',get(h1,'XTick'))
set(h5,'XTick',get(h1,'XTick'))
set(h6,'XTick',get(h1,'XTick'))
datetick(h2,'x',13,'keeplimits','keepticks')
datetick(h3,'x',13,'keeplimits','keepticks')
datetick(h4,'x',13,'keeplimits','keepticks')
datetick(h5,'x',13,'keeplimits','keepticks')
datetick(h6,'x',13,'keeplimits','keepticks')
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
set(cbh1,'Position',[cbsize1(1:2).*[1,1.02] .043716 .10319]);

cbsize2=get(cbh2,'Position');
set(cbh2,'Position',[cbsize2(1:2).*[1,1.02] .043716 .10319]);

cbsize3=get(cbh3,'Position');
set(cbh3,'Position',[cbsize3(1:2).*[1,1.02] .043716 .10319]);

cbsize4=get(cbh4,'Position');
set(cbh4,'Position',[cbsize4(1:2).*[1,1.02] .043716 .10319]);

cbsize5=get(cbh5,'Position');
set(cbh5,'Position',[cbsize5(1:2).*[1,1.02] .043716 .10319]);

set(h1,'fontsize',8)
set(cbh1,'fontsize',8)
set(h2,'fontsize',8)
set(cbh2,'fontsize',8)
set(h3,'fontsize',8)
set(cbh3,'fontsize',8)
set(h4,'fontsize',8)
set(cbh4,'fontsize',8)
set(h5,'fontsize',8)
set(cbh5,'fontsize',8)
set(h6,'fontsize',8)

if strcmp(p.Results.neplot,'log')
    set(cbh1,'YTick',5:13)
    set(cbh1,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})
    set(cbh2,'YTick',5:13)
    set(cbh2,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'}) ...

    set(cbh3,'YTick',5:13)
    set(cbh3,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'}) ...

end
set(cbh4,'YTick',5:13)
set(cbh4,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})


set(h1,'layer','top')
set(cbh1,'LineWidth',1,'TickDir','both')
set(h2,'layer','top')
set(cbh2,'LineWidth',1,'TickDir','both')
set(h3,'layer','top')
set(cbh3,'LineWidth',1,'TickDir','both')
set(h4,'layer','top')
set(cbh4,'LineWidth',1,'TickDir','both')
set(h5,'layer','top')
set(cbh5,'LineWidth',1,'TickDir','both')
set(h6,'layer','top')


tstr1 = datestr(datenum(datetime((ElSpecGlowOut.te(1)),'ConvertFrom', ...
                                 'posixtime')),29);
tstr2 = datestr(datenum(datetime((ElSpecGlowOut.te(end)),'ConvertFrom', ...
                                 'posixtime')),29);
if strcmp(tstr1,tstr2)
    tstr = tstr1;
else
    tstr = [tstr1,' -- ',tstr2];
end
title(h1,tstr)
drawnow

% adjust ylabel positions
posyh1 = get(yh1,'position');
posyh2 = get(yh2,'position');
posyh3 = get(yh3,'position');
posyh4 = get(yh4,'position');
posyh5 = get(yh5,'position');
posyh6 = get(yh6,'position');

set(yh1,'position',[posyh1(1) posyh1(2:3)]);
set(yh2,'position',[posyh1(1) posyh2(2:3)]);
set(yh3,'position',[posyh1(1) posyh3(2:3)]);
set(yh4,'position',[posyh1(1) posyh4(2:3)]);
set(yh5,'position',[posyh1(1) posyh5(2:3)]);
set(yh6,'position',[posyh1(1) posyh6(2:3)]);

posycbh1 = get(ycbh1,'position');
posycbh2 = get(ycbh2,'position');
posycbh3 = get(ycbh3,'position');
posycbh4 = get(ycbh4,'position');
posycbh5 = get(ycbh5,'position');

set(ycbh1,'position',[posycbh1(1) posycbh1(2:3)]);
set(ycbh2,'position',[posycbh1(1) posycbh2(2:3)]);
set(ycbh3,'position',[posycbh1(1) posycbh3(2:3)]);
set(ycbh4,'position',[posycbh1(1) posycbh4(2:3)]);
set(ycbh5,'position',[posycbh1(1) posycbh5(2:3)]);

h=findall(gcf,'tag','legend'); % find the handle to the legend
set(h,'location','none'); % Use the SET command to set the "Location" property

drawnow

end