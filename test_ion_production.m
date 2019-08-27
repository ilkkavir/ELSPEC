%
% test the two ion production calculations
%
%
%
% comments from tests:
%
% 2016-08-26:
%    there is a significant difference between the results of the
%    Sergienko & Ivanov model and those of the Fang et al. model.
%    One obvious reason for the difference is the Albedo flux
%    included the Sergienko model. Due to the albedo flux, about
%    half of the incident electron energy is not spent on ionizing
%    the neutrals, but is reflected back. By means of integrating
%    the ion production profiles with respect to height, it is easy
%    to see that the Fang et al. profiles use almost all
%    precipitating energy to ion pair production. Is this something
%    that is missing from the Fang model?!?!?
%
%
%
%
%
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later


H = 1000*(90:160);

[T rho] = atmosnrlmsise00(H, 70, 20, 2016, 20, 20*3600);

E=logspace(log10(100),log10(2e4),300);

nN2 = rho(:,3);
nO2 = rho(:,4);
nO = rho(:,2);
Tn = T(:,2);


[qS,EcS,dES] = ion_production_Sergienko1993(E,H,nN2,nO2,nO,Tn);
[qF,EcF,dEF] = ion_production_Fang2010(E,H,nN2,nO2,nO,Tn);

figure('Position', [100, 100, 500, 750]);
ax(1)=subplot(3,1,1);
pcolor(EcS/1000,H/1000,qS);shading flat;caxis([-1 1]* max([max(max(qS)) ...
                    max(max(qF))]))
title('Sergienko')
xlabel('Energy [keV]')
ylabel('Height [km]')

ax(2)=subplot(3,1,2);
pcolor(EcF/1000,H/1000,qF);shading flat;caxis([-1 1]* max([max(max(qS)) ...
                    max(max(qF))]));
title('Fang')
xlabel('Energy [keV]')
ylabel('Height [km]')

ax(3)=subplot(3,1,3);
pcolor(EcF/1000,H/1000,(qF-qS));shading flat;caxis([-1 1]* max([max(max(qS)) ...
                    max(max(qF))]));
title('Difference')
xlabel('Energy [keV]')
ylabel('Height [km]')

h=colorbar;

set(h, 'Position', [.8314 .11 .0581 .8150])
for i=1:3
    pos=get(ax(i), 'Position');
    set(ax(i), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
end

ylabel(h,'Ion production per precipitating electron [m^{-3}]')