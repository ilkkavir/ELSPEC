function photoprod = photoproductionGlow(ts,h,loc)
%
% photoionization rates with the glow model, modified from glow.m in the glowaurora package
%
%% GLOW model from Matlab.
% https://www.scivision.co/matlab-python-user-module-import/
assert(~verLessThan('matlab', '9.5'), 'Matlab >= R2018b required')

% change to the glow install directory to have the IRI files available
curpath = pwd;
%cd(fullfile(string(py.str(py.glowaurora.glowpath)),'glowaurora'));
cd(fileparts(which('glow.m')));

nt = length(ts);
nh = length(h);

photoprod = zeros(nh,nt);

for it=1:nt
    fprintf('\r %i / %i',it ,length(ts))
    
    timetmp = datetime(ts(it),'convertfrom','posixtime');
    pydate = py.datetime.datetime(int32(timetmp.Year),int32(timetmp.Month),int32(timetmp.Day),int32(timetmp.Hour),int32(timetmp.Minute));
    
    params = py.dict(pyargs('flux', 0, 'E0', 0, 'glat', loc(1), 'glon', loc(2), 't0', pydate));
    
    G = py.glowaurora.runglowaurora(params);
    
    photoion = xarray2mat(G{'photIon'});
    z_km = xarrayind2vector(G,'z_km');

    photoprod(:,it) = interp1(z_km,sum(photoion(:,1:2),2),h,'linear')*1e6;
    

end

% back to the original working directory
cd(curpath)

end


function M = xarray2mat(V)
M = double(py.numpy.asfortranarray(V));
end

function I = xarrayind2vector(V,key)
    
I = cell2mat(cell(V.indexes{key}.values.tolist)); 
    
end
