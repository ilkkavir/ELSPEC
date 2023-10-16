function photoprod = photoproductionGlow(ts,h,loc,f107,f107a,f107p,ap);
        %function photoprod = photoproductionGlow(ts,h,loc)
%
% photoionization rates with the glow model using glow.m in the ncar-glow package
%
% Using 177 energy bins in glow. This corresponds to maximum energy of 2.50e4 eV (see egrid.f90), which is larger than
% the photon energy of the shortest wavelength included (0.05 nm ~ 2.48e4 eV).
% 
%
%
% IV 2023

    nt = length(ts);
    nh = length(h);
    photoprod = NaN(nh,nt);



    % first calculate with 10-min time steps.
    ts2 = ts(1):600:ts(end);
    nt2 = length(ts2);
    photoprod2 = NaN(nh,nt2);
    f107a2 = interp1(ts,f107a,ts2);
    f1072 = interp1(ts,f107,ts2);
    f107p2 = interp1(ts,f107p,ts2);
    ap2 = interp1(ts,ap,ts2);
    for it2=1:nt2
        fprintf('\r %i / %i',it2 ,nt2)
        % if the glow run fails we skip it and interpolate later
        try
            glowout = ncarglow.glow(datetime(ts2(it2),'convertfrom','posixtime') , loc(1) , loc(2) , f107a2(it2) , f1072(it2) , f107p2(it2) , ap2(it2) , 0 , 1e2 , 177);
            
            photoprod2(:,it2) = interp1(glowout.altkm,glowout.ionrate,h,'linear')*1e6;
        catch
            ;
        end
    end
    
    % fill possible NaN values
    photoprod2 = fillmissing(photoprod2 , 'linear', 2 );
    photoprod2 = fillmissing(photoprod2 , 'linear', 1 );

    % interpolate to the final time resolution
    photoprod = interp2(ts2,h,photoprod2,ts,h,'spline');
    
end

    

% %% GLOW model from Matlab.
% % https://www.scivision.co/matlab-python-user-module-import/
% assert(~verLessThan('matlab', '9.5'), 'Matlab >= R2018b required')

% % change to the glow install directory to have the IRI files available
% curpath = pwd;
% %cd(fullfile(string(py.str(py.glowaurora.glowpath)),'glowaurora'));
% cd(fileparts(which('glow.m')));

% nt = length(ts);
% nh = length(h);

% photoprod = zeros(nh,nt);

% for it=1:nt
%     fprintf('\r %i / %i',it ,length(ts))
    
%     timetmp = datetime(ts(it),'convertfrom','posixtime');
%     pydate = py.datetime.datetime(int32(timetmp.Year),int32(timetmp.Month),int32(timetmp.Day),int32(timetmp.Hour),int32(timetmp.Minute));
    
%     params = py.dict(pyargs('flux', 0, 'E0', 0, 'glat', loc(1), 'glon', loc(2), 't0', pydate));
    
%     G = py.glowaurora.runglowaurora(params);
    
%     photoion = xarray2mat(G{'photIon'});
%     z_km = xarrayind2vector(G,'z_km');

%     photoprod(:,it) = interp1(z_km,sum(photoion(:,1:2),2),h,'linear')*1e6;
    

% end

% % back to the original working directory
% cd(curpath)

% end


% function M = xarray2mat(V)
% M = double(py.numpy.asfortranarray(V));
% end

% function I = xarrayind2vector(V,key)
    
% I = cell2mat(cell(V.indexes{key}.values.tolist)); 
    
% end
