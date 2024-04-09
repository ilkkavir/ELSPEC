function ElSpecGlowOut = GLOWemissions(ElSpecOut,saveresults)
%
% ElSpecGlowOut = GLOWemissions(ElSpecOut)
%
% Calculate emission profiles using spectra from ELSPEC and the GLOW model.
%
% The program glowELSPEC must be on the search path under the name glowelspec
%
%
% IV 2020
%

    if nargin==1
        saveresults = false;
    end

    Pmin = 1e-5;
    Epeakmin = 8e2;
    
% temporary file names
    glowin = fullfile(tempdir,'glowinput.dat');
    glowout = fullfile(tempdir,'glowoutput.dat');
    
% timestamps as matlab datetime
    times = datetime( (ElSpecOut.ts + ElSpecOut.te)/2 ,'convertfrom','posixtime');
    yyyy = times.Year;
    ddd = floor(datenum(times)-datenum(dateshift(times,'start','year'))) + 1;
    utsec = seconds(times-dateshift(times,'start','day'));


    % write the glow inputs in file, one row per time step
    nt = length(ElSpecOut.ts);
    fidin = fopen(glowin,'w');
    for it = 1:nt
        % coefficients of the polynomial model, assume 10 coefficients, pad with zeros as necessary
        polycoefs = zeros(10,1);
        if ~isnan(ElSpecOut.best_order(it))
            for icoef = 1:(ElSpecOut.best_order(it)+1)
                polycoefs(icoef) = ElSpecOut.polycoefs(ElSpecOut.best_order(it),icoef,it);
            end
        end
        % write to file
        fprintf(fidin,'%04d%03d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',yyyy(it),ddd(it),utsec(it),ElSpecOut.loc(1),ElSpecOut.loc(2),ElSpecOut.f107a(it),ElSpecOut.f107(it),ElSpecOut.f107p(it),ElSpecOut.ap(it),polycoefs(1),polycoefs(2),polycoefs(3),polycoefs(4),polycoefs(5),polycoefs(6),polycoefs(7),polycoefs(8),polycoefs(9),polycoefs(10));
    end
    fclose(fidin);

    system(['glowelspec <' glowin '> ' glowout]);

    delete(glowin);

    nh = 102; % this is set in glowELSPEC.f90
    
    ElSpecGlowOut = ElSpecOut;
    ElSpecGlowOut.GLOWh = NaN(nh,1);
    ElSpecGlowOut.GLOWne = NaN(nh,nt);
    ElSpecGlowOut.GLOWe5577 = NaN(nh,nt);
    ElSpecGlowOut.GLOWe4278 = NaN(nh,nt);
    ElSpecGlowOut.GLOWtime = times;

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
                ElSpecGlowOut.GLOWh(ih) = str2num(strparts{2});
            end
            ElSpecGlowOut.GLOWne(ih,it) = str2num(strparts{8})*1e6;
            ElSpecGlowOut.GLOWe5577(ih,it) = str2num(strparts{19});
            ElSpecGlowOut.GLOWe4278(ih,it) = str2num(strparts{17});
        end
        ih = ih + 1;
    end

    fclose(fido);

    delete(glowout)

    [ElSpecGlowOut.GLOWe5577max i5577max] = max(ElSpecGlowOut.GLOWe5577);
    ElSpecGlowOut.GLOWh5577max = ElSpecGlowOut.GLOWh(i5577max);
    [ElSpecGlowOut.GLOWe4278max i4278max] = max(ElSpecGlowOut.GLOWe4278);
    ElSpecGlowOut.GLOWh4278max = ElSpecGlowOut.GLOWh(i4278max);
    [ElSpecGlowOut.nemax inemax] = max(ElSpecGlowOut.ne);
    ElSpecGlowOut.hnemax = ElSpecGlowOut.h(inemax);
    glownetmp = ElSpecGlowOut.GLOWne;
    glownetmp(ElSpecGlowOut.GLOWh>150,:) = 0;
    [ElSpecGlowOut.GLOWnemax GLOWinemax] = max(glownetmp);
    ElSpecGlowOut.GLOWhnemax = ElSpecGlowOut.GLOWh(GLOWinemax);
    dh = diff(ElSpecGlowOut.GLOWh);
    dh = [dh(1) ; dh];
    ElSpecGlowOut.GLOWh5577mean = sum(ElSpecGlowOut.GLOWh.*dh.*ElSpecGlowOut.GLOWe5577,'omitnan')./sum(dh.*ElSpecGlowOut.GLOWe5577,'omitnan');
    ElSpecGlowOut.GLOWh4278mean = sum(ElSpecGlowOut.GLOWh.*dh.*ElSpecGlowOut.GLOWe4278,'omitnan')./sum(dh.*ElSpecGlowOut.GLOWe4278,'omitnan');

    % peak energy of differential number flux
    [Iepeak indIepeak] = max(ElSpecGlowOut.Ie);
    Ecpeak = ElSpecGlowOut.Ec(indIepeak);

    %    rminds = ElSpecGlowOut.Pe < Pmin | Ecpeak < Epeakmin;
    rminds = ElSpecGlowOut.Pe < Pmin ;%| Ecpeak < Epeakmin;
    ElSpecGlowOut.GLOWh5577max(rminds) = NaN;
    ElSpecGlowOut.GLOWe5577max(rminds) = NaN;
    ElSpecGlowOut.GLOWh5577mean(rminds) = NaN;
    ElSpecGlowOut.GLOWe5577(:,rminds) = NaN;
    ElSpecGlowOut.GLOWh4278max(rminds) = NaN;
    ElSpecGlowOut.GLOWe4278max(rminds) = NaN;
    ElSpecGlowOut.GLOWh4278mean(rminds) = NaN;
    ElSpecGlowOut.GLOWe4278(:,rminds) = NaN;

    if saveresults
        outfilename = ['GLOW_ElSpec_',...
                       datestr(datetime(round(ElSpecOut.ts(1)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),...
                       '-',...
                       datestr(datetime(round(ElSpecOut.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),...
                       '_',...
                       ElSpecOut.experiment,'_',ElSpecOut.radar,'_',ElSpecOut.ionomodel,'_',ElSpecOut.recombmodel,'_',ElSpecOut.integtype,...
                       '_',num2str(ElSpecOut.ninteg),'_',num2str(ElSpecOut.nstep),'_',ElSpecOut.tres,'_',...
                       datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];
        
        try
            save(outfilename,'ElSpecGlowOut','-v7.3');
        catch
            ;
        end
    end

end