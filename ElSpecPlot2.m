function fignum = ElSpecPlot2(ElSpecOut,varargin)
%
% plot ion production rate, recombination time-scale, dNe/dT,
% instead of flux, FAC, power
%
%
%
%
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

defaultProdlim = [8 12];
checkProdlim = @(x) (isnumeric(x) & (length(x)==2));

defaultRecomblim = [0 2];
checkRecomblim = @(x) (isnumeric(x) & (length(x)==2));

defaultDnelim = [-1 1]*1e10;
checkDnelim = @(x) (isnumeric(x) & (length(x)==2));

defaultChisqrlim = [0 5];
checkChisqrlim = @(x) (isnumeric(x) & (length(x)==2));

defaultBtime = datetime((ElSpecOut.te(1) + ElSpecOut.ts(1))/2,'ConvertFrom','posixtime');
checkBtime = @(x) (isnumeric(x) & (length(x)==6));

defaultEtime = datetime((ElSpecOut.te(end) + ElSpecOut.ts(end))/2,'ConvertFrom','posixtime');
checkEtime = @(x) (isnumeric(x) & (length(x)==6));

defaultNeplot = 'log';
validNeplot = {'log','linear'};
checkNeplot = @(x) any(validatestring(x,validNeplot));

addRequired(p,'ElSpecOut',@isstruct);
addOptional(p,'fignum',defaultFignum,@isnumeric);
addParameter(p,'hlim',defaultHlim,checkHlim)
addParameter(p,'nelim',defaultNelim,checkNelim)
addParameter(p,'prodlim',defaultProdlim,checkProdlim)
addParameter(p,'recomblim',defaultRecomblim,checkRecomblim)
addParameter(p,'dnelim',defaultDnelim,checkDnelim)
addParameter(p,'chisqrlim',defaultChisqrlim,checkChisqrlim)
addParameter(p,'btime',defaultBtime,checkBtime)
addParameter(p,'etime',defaultEtime,checkEtime)
addParameter(p,'neplot',defaultNeplot,checkNeplot)
parse(p,ElSpecOut,varargin{:})

% recombination timescale
recombscale = 1./(ElSpecOut.alpha.*ElSpecOut.ne);

% electron density time-derivative
dne = diff(ElSpecOut.ne')';
dne(:,end+1) = 0;
for it = 1:length(ElSpecOut.dt)
    dne(:,it) = dne(:,it) ./ ElSpecOut.dt(it);
end

% centres of time-bins as matlab datenum
tt = datenum(datetime((ElSpecOut.te + ElSpecOut.ts)/2,'ConvertFrom','posixtime'));
% start points of time-bins for the pcolor plots
ts = datenum(datetime(ElSpecOut.ts,'ConvertFrom','posixtime'));
% heighs
hh = ElSpecOut.h;
% power profiles with negatives removed
pp = ElSpecOut.pp;
pp(pp<1) = 1;

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
set(gcf,'Position',figpos,'PaperPositionMode','Auto');
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
ylabel('Height [km]')
cbh1=colorbar;
ylabel(cbh1,{'N_e [m^{-3}] (meas.)'})

h2=subplot(6,1,2);
if strcmp(p.Results.neplot,'log')
    pcolor(ts,hh,log10(ElSpecOut.ne)),shading flat
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

h3=subplot(6,1,3);
pcolor(ts,hh,log10(ElSpecOut.q)),shading flat
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(p.Results.prodlim)
ylabel('Height [km]')
cbh3=colorbar;
ylabel(cbh3,'Production [m^{-3}s^{-1}]')
%ylim([.5 20])

h4=subplot(6,1,4);
pcolor(ts,hh,log10(recombscale)),shading flat
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
%caxis(p.Results.ielim)
caxis(p.Results.recomblim)
ylabel('Height [km]')
cbh4=colorbar;
ylabel(cbh4,'Rec. time [s]')


h5=subplot(6,1,5);
pcolor(ts,hh,dne),shading flat
ylim(p.Results.hlim)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
caxis(p.Results.dnelim)
ylabel('Height [km]')
cbh5=colorbar;
ylabel(cbh5,'DN_e [m^{-3}s^{-1}]')

h6=subplot(6,1,6);
plot(tt,ElSpecOut.chisqr)
xlim(datenum([datetime(p.Results.btime) datetime(p.Results.etime)]))
ylim(p.Results.chisqrlim);
ylabel('\chi^2')
xlabel('Time [UTC]');


drawnow % this is important to update everything before setting the sizes...
pos1=get(h1,'Position');
pos2=get(h2,'Position');
pos3=get(h3,'Position');
pos4=get(h4,'Position');
pos5=get(h5,'Position');
pos6=get(h6,'Position');
set(h1, 'Position', [pos1(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
set(h2, 'Position', [pos2(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
set(h3, 'Position', [pos3(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
set(h4, 'Position', [pos4(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
set(h5, 'Position', [pos5(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
set(h6, 'Position', [pos6(1:2)-[0,.05] [.85 1.77].*pos4(3:4)]);
linkaxes([h1 h2 h3 h4 h5 h6],'x')
linkaxes([h1 h2 h3 h4],'y')

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
set(h3,'TickDir','both')
set(h3,'YTick',hticks)
set(h4,'TickDir','both')
set(h4,'YTick',hticks)
set(h5,'TickDir','both')
set(h5,'YTick',hticks)
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
set(cbh1,'Position',cbsize1.*[1,1.02,1,.8]);

cbsize2=get(cbh2,'Position');
set(cbh2,'Position',cbsize2.*[1,1.02,1,.8]);

cbsize3=get(cbh3,'Position');
set(cbh3,'Position',cbsize3.*[1,1.02,1,.8]);

cbsize4=get(cbh4,'Position');
set(cbh4,'Position',cbsize4.*[1,1.02,1,.8]);

cbsize5=get(cbh5,'Position');
set(cbh5,'Position',cbsize5.*[1,1.02,1,.8]);

set(h1,'fontsize',12)
set(cbh1,'fontsize',12)
set(h2,'fontsize',12)
set(cbh2,'fontsize',12)
set(h3,'fontsize',12)
set(cbh3,'fontsize',12)
set(h4,'fontsize',12)
set(cbh4,'fontsize',12)
set(h5,'fontsize',12)
set(cbh5,'fontsize',12)
set(h6,'fontsize',12)

if strcmp(p.Results.neplot,'log')
    set(cbh1,'YTick',5:13)
    set(cbh1,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'})
    set(cbh2,'YTick',5:13)
    set(cbh2,'YTickLabel',{'10^{5}';'10^{6}';'10^{7}';'10^{8}'; ...
                        '10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}'}) ...

end

set(cbh3,'YTick',6:15)
set(cbh3,'YTickLabel',{'10^{6}';'10^{7}';'10^{8}';'10^{9}';'10^{10}';'10^{11}';'10^{12}';'10^{13}';'10^{14}';'10^{15}'})
set(cbh4,'YTick',-1:3)
set(cbh4,'YTickLabel',{'0.1';'1';'10';'100';'1000'})
set(cbh5,'YTick',[-1 0 1]*1e10)
set(cbh5,'YTickLabel',{'-10^{10}';'0';'10^{10}'})


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