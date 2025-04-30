function createModelParameterFile2(varargin)
%
%
% Tabulate the necessary IRI parameters in Tromso and Svalbard
% We will need O2+, NO+, and Te, but also O+ and Ti are tabulated
% because they may be very useful outside ELSPEC
%
% Converted to MATLAB from the original R version, 2025. 
%
% The function uses IRI model from GUISDAP. The easiest way to set the matlab path correctly is to
% start matlab from the command line with the command 'guisdap'. 
%
%
% IV 2017, 2022, 2025
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi>
% This is free software, licensed under GNU GPL version 2 or later

    
% parse the input
    p = inputParser;
    
    % Start end years
    defaultStartYear = 2000;
    defaultEndYear = datetime('now').Year;
    checkYear = @(x) (isnumeric(x) & length(x)==1 & x >= 1980 );

    % height limits
    defaultHmin = 80;
    defaultHmax = 150;
    checkH = @(x) (isnumeric(x) & length(x)==1 & x > 0 & x < 1000);

    % height resolution
    defaultHres = 2;
    checkHres = @(x) (isnumeric(x) & length(x)==1 & x > 0 & x < 10);

    % output file name
    defaultFname = 'ElSpecModelParameters';
    checkFname = @(x) (isstr(x));
    
    addParameter(p,'startyear',defaultStartYear,checkYear);
    addParameter(p,'endyear',defaultEndYear,checkYear);
    addParameter(p,'hmin',defaultHmin,checkH);
    addParameter(p,'hmax',defaultHmax,checkH);
    addParameter(p,'hres',defaultHres,checkHres);
    addParameter(p,'fname',defaultFname,checkFname);
    parse(p,varargin{:})


    years = p.Results.startyear:p.Results.endyear;
    nyear = length(years);

    h = p.Results.hmin:p.Results.hres:p.Results.hmax;
    nh = length(h);

    ntot = nyear*12*31*24;

    ns = 0;
    for iyear = 1:nyear
        
        iriparTRO = NaN(12,31,24,nh,5);
        iriparESR = NaN(12,31,24,nh,5);
        
        for month = 1:12
            for day = 1:31
                for hour = 0:23
                    ns = ns+1;
                    tcur = datetime(years(iyear),month,day,hour,0,0);
                    tsec = seconds(tcur - dateshift(tcur,'start','year'));
                    if (years(iyear)>1980)
                        iriparTRO(month,day,hour+1,:,:) = iri([5,8,9,3,4],[tsec,years(iyear)],[69.5864,19.2272],[min(h) max(h) p.Results.hres]);
                    end
                    if years(iyear)>1995
                        iriparESR(month,day,hour+1,:,:) = iri([5,8,9,3,4],[tsec,years(iyear)],[78.15,16.02],[min(h) max(h) p.Results.hres]);                        
                    end
                    fprintf("                 \r")
                    fprintf("%6.0f %2.0f %2.0f %2.0f %3.0f%%",years(iyear),month,day,hour,ns/ntot*100)
                end
            end
        end
        if years(iyear)>1980
            save([p.Results.fname,'TRO',num2str(years(iyear)),'.mat'],'iriparTRO');
        end
        if years(iyear)>1995
            save([p.Results.fname,'ESR',num2str(years(iyear)),'.mat'],'iriparESR');
        end
    end
end
