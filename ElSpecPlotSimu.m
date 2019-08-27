function fignum = ElSpecPlotSimu(ElSpecOut,simudata,varargin)
%
% Plot simulation result. This is ElSpecPlot with an additional
% panel for the known spectrum used in the simulation. The known
% true power and true FAC are plotted on top of the fit results.
%
% IV 2017
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later



%fignum = -1;
%if nargin > 1
%    fignum = varargin{1};
%end

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

defaultElim = [1 50];%[min(ElSpecOut.Ec) max(ElSpecOut.Ec)]./1000;
checkElim = @(x) (isnumeric(x) & (length(x)==2));

defaultFAClim = [0 10];
checkFAClim = @(x) (isnumeric(x) & (length(x)==2));

defaultPlim = [0 30];
checkPlim = @(x) (isnumeric(x) & (length(x)==2));

defaultChisqrlim = [0 5];
checkChisqrlim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecOut.te(1) + ElSpecOut.ts(1))/2 - ElSpecOut.ts(1),'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecOut.te(end) + ElSpecOut.ts(end))/2 ...
                        - ElSpecOut.ts(1),'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultFluxtype = 'energy';
validFluxtype = {'number','energy'};
checkFluxtype = @(x) any(validatestring(x,validFluxtype));

defaultNeplot = 'log';
validNeplot = {'log','linear'};
checkNeplot = @(x) any(validatestring(x,validNeplot));


if isfield(ElSpecOut,'Emin')
    defaultEmin = ElSpecOut.Emin;
else
    defaultEmin = 1e3;
end
checkEmin = @(x) (isnumeric(x) & length(x)==1 );

defaultLabelStep = 0;
checkLabelStep = @(x) (isnumeric(x) & (length(x)==1));

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
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'fluxtype',defaultFluxtype,checkFluxtype)
addParameter(p,'neplot',defaultNeplot,checkNeplot)
addParameter(p,'emin',defaultEmin,checkEmin)
addParameter(p,'labelstep',defaultLabelStep,checkLabelStep)
parse(p,ElSpecOut,varargin{:})


% centres of time-bins as matlab datenum, zero-time is at beginning
% of the simulation
tt = datenum(datetime((ElSpecOut.te + ElSpecOut.ts)/2-ElSpecOut.ts(1),'ConvertFrom','posixtime'));
% start points of time-bins for the pcolor plots
ts = datenum(datetime(ElSpecOut.ts-ElSpecOut.ts(1),'ConvertFrom','posixtime'));

% heighs
hh = ElSpecOut.h;
% energies
EE = ElSpecOut.Ec/1000;
EEsimu = simudata.Ec/1000;
% power profiles with negatives removed
pp = ElSpecOut.pp;
pp(pp<1) = 1;

% update FAC and power using the Emin input
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
Eind_fac_simu = simudata.Ec >= p.Results.emin;
for tind=1:length(simudata.ts)
    simudata.FAC(tind) = sum(simudata.Ie(Eind_fac_simu,tind).* ...
                              simudata.dE(Eind_fac_simu)')*1.60217662e-19;
    % Power carried by the precipitating electrons
    simudata.Pe(tind) = sum(simudata.Ie(Eind_fac_simu,tind).*simudata.dE(Eind_fac_simu)'.*simudata.Ec(Eind_fac_simu)')*1.60217662e-19;
    EdE = simudata.dE(Eind_fac_simu)'.*simudata.Ec(Eind_fac_simu)';

end

% number vs energy flux
switch lower(p.Results.fluxtype)

  case 'number'
    Ie = ElSpecOut.Ie;
    Ielabel={'[(eVm^2s)^{-1}]'};
    ielim = p.Results.ielim;
    IeSimu = simudata.Ie
  case 'energy'
    Ie = ElSpecOut.Ie;
    IeSimu = simudata.Ie;
    for iiE = 1:length(ElSpecOut.Ec)
        Ie(iiE,:) = Ie(iiE,:).*ElSpecOut.Ec(iiE);
    end
    for iiE = 1:length(simudata.Ec)
        IeSimu(iiE,:) = IeSimu(iiE,:).*simudata.Ec(iiE);
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
    %    figpos = [0 0 21/29.7*690 690];
    figpos = [0 0 21/29.7*690 690];
end
%figpos(3:4) = [21/29.7*690 690];
figpos(3:4) = [21/29.7*690 690];
set(gcf,'Position',figpos,'PaperPositionMode','Auto');
colormap(jet);

h1=subplot(7,1,1);
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
ylabel('Height [km]')
cbh1=colorbar;
ylabel(cbh1,{'N_e [m^{-3}] (simul.)'})

h2=subplot(7,1,2);
if strcmp(p.Results.neplot,'log')
    ppOut = ElSpecOut.ne;
    ppOut(ppOut<1) = 1;
    pcolor(ts,hh,log10(ppOut)),shading flat
    caxis(p.Results.nelim)
elseif strcmp(p.Results.neplot,'linear')
    pcolor(ts,hh,ElSpecOut.ne),shading flat
    caxis(10.^(p.Results.nelim))
else
    error(['Unknown neplot ',neplot])
end
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.hlim)
cbh2=colorbar;
ylabel('Height [km]')
ylabel(cbh2,{'N_e [m^{-3}] (model)'})

h3=subplot(7,1,3);
pcolor(ts,EEsimu,log10(IeSimu)),shading flat
set(h3,'yscale','log')
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(ielim)
ylabel('Energy [keV]')
cbh3=colorbar;
ylabel(cbh3,Ielabel)
%ylim([.5 20])

h4=subplot(7,1,4);
pcolor(ts,EE,log10(Ie)),shading flat
set(h4,'yscale','log')
ylim(p.Results.elim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(ielim)
ylabel('Energy [keV]')
cbh4=colorbar;
ylabel(cbh4,Ielabel)
%ylim([.5 20])


h5=subplot(7,1,5);
plot(tt,ElSpecOut.FAC*NaN)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.faclim)
hold on
for tind=1:length(ElSpecOut.ts)
    plot([tt(tind) tt(tind)] , [-1 1]*ElSpecOut.FACstd(tind)*1e6+ElSpecOut.FAC(tind)*1e6,'r')
end
plot(tt,ElSpecOut.FAC*1e6,'k-','LineWidth',1)
plot(tt,simudata.FAC*1e6,'b-','LineWidth',1)
hold off
ylabel('FAC [\muAm^{-2}]')
grid on


h6=subplot(7,1,6);
plot(tt,ElSpecOut.Pe*NaN)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.plim)
hold on
for tind=1:length(ElSpecOut.ts)
    plot([tt(tind) tt(tind)] , [-1 1]*ElSpecOut.PeStd(tind)*1e3+ElSpecOut.Pe(tind)*1e3,'r')
end
plot(tt,ElSpecOut.Pe*1000,'k-','LineWidth',1)
plot(tt,simudata.Pe*1e3,'b-','LineWidth',1)
hold off
ylabel('Power [mWm^{-2}]')
grid on

h7=subplot(7,1,7);
plot(tt,ElSpecOut.chisqr,'k-','LineWidth',1)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.chisqrlim);
ylabel('\chi^2')
xlabel('Time');
grid on


drawnow % this is important to update everything before setting the sizes...
pos1=get(h1,'Position');
pos2=get(h2,'Position');
pos3=get(h3,'Position');
pos4=get(h4,'Position');
pos5=get(h5,'Position');
pos6=get(h6,'Position');
pos7=get(h7,'Position');
set(h1, 'Position', [pos1(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h2, 'Position', [pos2(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h3, 'Position', [pos3(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h4, 'Position', [pos4(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h5, 'Position', [pos5(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h6, 'Position', [pos6(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
set(h7, 'Position', [pos7(1:2)-[0,.05] [.85 1.98].*pos5(3:4)]);
linkaxes([h1 h2 h3 h4 h5 h6 h7],'x')
linkaxes([h1 h2],'y')
linkaxes([h3 h4],'y')

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
set(h4,'YTick',eticks(eticks<p.Results.elim(2)))
set(h4,'TickDir','both')
set(h5,'TickDir','both')
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
set(h5,'YTick',facticks)
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
set(h6,'TickDir','both')
chistep=1;
if range(p.Results.chisqrlim) > 5
    chistep=2;
end
if range(p.Results.chisqrlim) > 10
    chistep=4;
end
chiticks = ceil(p.Results.chisqrlim(1)):chistep:floor(p.Results.chisqrlim(end)-.1);
set(h7,'TickDir','both')
set(h7,'Ytick',chiticks)
set(h2,'XTick',get(h1,'XTick'))
set(h3,'XTick',get(h1,'XTick'))
set(h4,'XTick',get(h1,'XTick'))
set(h5,'XTick',get(h1,'XTick'))
set(h6,'XTick',get(h1,'XTick'))
set(h7,'XTick',get(h1,'XTick'))
datetick(h2,'x',13,'keeplimits','keepticks')
datetick(h3,'x',13,'keeplimits','keepticks')
datetick(h4,'x',13,'keeplimits','keepticks')
datetick(h5,'x',13,'keeplimits','keepticks')
datetick(h6,'x',13,'keeplimits','keepticks')
datetick(h7,'x',13,'keeplimits','keepticks')
set(h1,'XTickLabel','')
set(h2,'XTickLabel','')
set(h3,'XTickLabel','')
set(h4,'XTickLabel','')
set(h5,'XTickLabel','')
set(h6,'XTickLabel','')
set(h1,'LineWidth',1)
set(h2,'LineWidth',1)
set(h3,'LineWidth',1)
set(h4,'LineWidth',1)
set(h5,'LineWidth',1)
set(h6,'LineWidth',1)
set(h7,'LineWidth',1)

cbsize1=get(cbh1,'Position');
set(cbh1,'Position',cbsize1.*[1,1.02,1,.8]);

cbsize2=get(cbh2,'Position');
set(cbh2,'Position',cbsize2.*[1,1.02,1,.8]);

cbsize3=get(cbh3,'Position');
set(cbh3,'Position',cbsize3.*[1,1.02,1,.8]);

cbsize4=get(cbh4,'Position');
set(cbh4,'Position',cbsize4.*[1,1.02,1,.8]);

set(h1,'fontsize',12)
set(cbh1,'fontsize',12)
set(h2,'fontsize',12)
set(cbh2,'fontsize',12)
set(h3,'fontsize',12)
set(cbh3,'fontsize',12)
set(h4,'fontsize',12)
set(cbh4,'fontsize',12)
set(h5,'fontsize',12)
set(h6,'fontsize',12)
set(h7,'fontsize',12)

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
set(h6,'layer','top')
set(h7,'layer','top')


tstr1 = datestr(datenum(datetime((ElSpecOut.te(1)),'ConvertFrom', ...
                                 'posixtime')),29);
tstr2 = datestr(datenum(datetime((ElSpecOut.te(end)),'ConvertFrom', ...
                                 'posixtime')),29);
if strcmp(tstr1,tstr2)
    tstr = tstr1;
else
    tstr = [tstr1,' -- ',tstr2];
end
title(h1,'')
drawnow

% optional labels A,B,C,... above the plot
if p.Results.labelstep > 0
    hlab = axes('Position',get(h1,'Position'),'box','off','xtick',[],'ytick',[],'color','none');
    set(hlab,'xaxisLocation','top')
    set(hlab,'XLim',get(h1,'XLim'))
    xlims = get(h1,'XLim');
    % not sure why we need to subtract 1...
    labticks = xlims(1)-1/24/3600 : p.Results.labelstep/3600/24 ...
        : xlims(2);
    Alphabet = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'};
    set(hlab,'XTick',labticks)
    set(hlab,'XTickLabel',[])
    set(hlab,'layer','bottom')
    set(hlab,'LineWidth',1,'fontsize',12,'TickDir','out')
    hlab2 = axes('Position',get(hlab,'Position'),'box','off','xtick',[],'ytick',[],'color','none','TickLength',[0,0]);
    set(hlab2,'XLim',get(hlab,'XLim'))
    set(hlab2,'xaxisLocation','top')
    set(hlab2,'XTick',labticks+(labticks(2)-labticks(1))/2)
    set(hlab2,'XTickLabel',Alphabet)
    set(hlab2,'LineWidth',1,'fontsize',12,'TickDir','out')
    drawnow
end