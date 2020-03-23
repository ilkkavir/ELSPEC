function GlowOut = peakNepeakEmissionTest( time , Q0 , E0)
%
% GlowOut = peakNepeakEmissionTest( time , Q0 , E0 )
%
% test the correlation between peak Ne and peak emission altitudes using different maxwellian fluxes as input
%
% INPUT:
%    time    time as matlab datetime
%    Q0      total energy flux (mW/m^2)
%    E0      a vector of characteristic energies (keV)
%
%
%
% The program glowbasic must be on the search path under the name glowbasic
%
%
% IV 2020
%


    nE = length(E0);
% indices and location
    f107a = 70;
    f107 = 70;
    f107p = 70;
    ap = 12;
    loc = [ 69.583  19.21  0.03];

    % temporary file names
    glowin = fullfile(tempdir,'glowinput.dat');
    glowout = fullfile(tempdir,'glowoutput.dat');
    
    % timestamp
    yyyy = time.Year;
    ddd = floor(datenum(time)-datenum(dateshift(time,'start','year'))) + 1;
    utsec = seconds(time-dateshift(time,'start','day'));


    % write the glow inputs in file
    fidin = fopen(glowin,'w');
    for k=1:nE
        fprintf(fidin,'%04d%03d %f %f %f %f %f %f %f %f %f\n',yyyy,ddd,utsec,loc(1),loc(2),f107a,f107,f107p,ap,Q0/1.60217662e-19/6.24150913e11/1e4/1e3,E0(k)*1e3);
    end
    fclose(fidin);

    system(['glowelspecmaxwell <' glowin '> ' glowout]);

    delete(glowin);

    nh = 102; % this is set in glowELSPEC.f90
    
    GlowOut = struct();
    GlowOut.GLOWh = NaN(nh,1);
    GlowOut.GLOWne = NaN(nh,nE);
    GlowOut.GLOWe5577 = NaN(nh,nE);
    GlowOut.GLOWtime = time;

    fido = fopen(glowout,'r');

    it = 0;
    ih = -10;
    while ~feof(fido)
        fline = fgetl(fido);
        if strcmp(fline(1:10),'##########')
            it = it + 1;
            ih = -2;
        end
        if strcmp(fline(1:10),'----------')
            ih = -10;
        end
        if ih > 0
            strparts = strsplit(fline);
            if it==1
                GlowOut.GLOWh(ih) = str2num(strparts{2});
            end
            GlowOut.GLOWne(ih,it) = str2num(strparts{8})*1e6;
            GlowOut.GLOWe5577(ih,it) = str2num(strparts{19});
        end
        ih = ih + 1;
    end

    fclose(fido);

    delete(glowout)

    [GlowOut.GLOWe5577max i5577max] = max(GlowOut.GLOWe5577);
    GlowOut.GLOWh5577max = GlowOut.GLOWh(i5577max);
    glownetmp = GlowOut.GLOWne;
    glownetmp(GlowOut.GLOWh>150,:) = 0;
    [GlowOut.GLOWnemax GLOWinemax] = max(glownetmp);
    GlowOut.GLOWhnemax = GlowOut.GLOWh(GLOWinemax);

    GlowOut.E0 = E0;

    f1 = figure;
    set(f1,'position',[0 0 800 800])
    h1=subplot(3,1,1);
    pcolor(GlowOut.E0,GlowOut.GLOWh,GlowOut.GLOWne)
    shading flat
    cb1=colorbar;
    ylim([80 150])
    xlabel('E0 [keV]')
    ylabel('Altitude [km]')
    ylabel(cb1,'Ne [m^{-3}]')
    title(['Auroral power = ' num2str(Q0) ' mW/m^2'])

    h2=subplot(3,1,2);
    pcolor(GlowOut.E0,GlowOut.GLOWh,GlowOut.GLOWe5577)
    shading flat
    cb1=colorbar;
    ylim([80 150])
    xlabel('E0 [keV]')
    ylabel('Altitude [km]')
    ylabel(cb1,'557.7 nm [cm^{-3}s^{-1}]')

    h3=subplot(3,1,3);
    plot(GlowOut.E0,GlowOut.GLOWh5577max,'linewidth',1.5)
    hold on
    plot(GlowOut.E0,GlowOut.GLOWhnemax,'linewidth',1.5)
    xlabel('E0 [keV]')
    ylabel('Altitude [km]')
    ylim([80 150])
    legend('peak emission','peak N_e')
    grid on

    pos1 = get(h1,'position');
    pos2 = get(h2,'position');
    pos3 = get(h3,'position');

    set(h3,'position',[pos3(1:2) pos1(3) pos3(4)])

   

end