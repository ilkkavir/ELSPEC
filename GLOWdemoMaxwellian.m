function GlowOut = GLOWdemoMaxwellian( time , Q0 , E0 )
%
% GlowOut = GLOWdemoMaxwellian( time , Q0 , E0)
%
% Calculate emission profiles using Maxwellian spectra from ELSPEC and the GLOW model.
%
% INPUT:
%    time    time as matlab datetime
%    Q0      total energy flux (mW/m^2)
%    E0      characteristic energy (keV)
%
%
%
% The program glowbasic must be on the search path under the name glowbasic
%
%
% IV 2020
%


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
    fprintf(fidin,'%04d%03d %f %f %f %f %f %f %f %f %f\n',yyyy,ddd,utsec,loc(1),loc(2),f107a,f107,f107p,ap,Q0/1.60217662e-19/6.24150913e11/1e4/1e3,E0*1e3);
    fclose(fidin);

    system(['glowelspecmaxwell <' glowin '> ' glowout]);

    delete(glowin);

    nh = 102; % this is set in glowELSPEC.f90
    
    GlowOut = struct();
    GlowOut.GLOWh = NaN(nh,1);
    GlowOut.GLOWne = NaN(nh,1);
    GlowOut.GLOWe5577 = NaN(nh,1);
    GlowOut.GLOWtime = time;

    fido = fopen(glowout,'r');

    ih = -10;
    while ~feof(fido)
        fline = fgetl(fido);
        if strcmp(fline(1:10),'##########')
            ih = -2;
        end
        if strcmp(fline(1:10),'----------')
            ih = -10;
        end
        if ih > 0
            strparts = strsplit(fline);
            GlowOut.GLOWh(ih) = str2num(strparts{2});
            GlowOut.GLOWne(ih) = str2num(strparts{8})*1e6;
            GlowOut.GLOWe5577(ih) = str2num(strparts{19});
        end
        ih = ih + 1;
    end

    fclose(fido);

    delete(glowout)

    [GlowOut.GLOWe5577max i5577max] = max(GlowOut.GLOWe5577);
    GlowOut.GLOWh5577max = GlowOut.GLOWh(i5577max);
    dh = diff(GlowOut.GLOWh);
    dh = [dh(1) ; dh];
    GlowOut.GLOWh5577mean = sum(GlowOut.GLOWh.*dh.*GlowOut.GLOWe5577,'omitnan')./sum(dh.*GlowOut.GLOWe5577,'omitnan');


end