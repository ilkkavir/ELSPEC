function modelout = calculateFlipchemComposition(ts,h,par,pp,loc,modelin)
%
% model = calculateFlipchemComposition(ts,h,par,loc)
%
% Calculate  NO+ and O2+ fractions using the flipchem model
%
% INPUT:
%   ts      integretion period start times [unixtime]
%   h       altitudes [km9
%   par     a parameter array as returned by readFitData
%   loc     location as returned by readFitData
%   modelin an array of model parameters as returned by readFitData
%
% OUTPUT:
%   modelout  an array of model parameters with ion compositions replaced with the flipchem values
%
%
% IV 2021

    modelout = modelin;
    for it = 1:length(ts)
        fprintf('\r %i / %i',it ,length(ts))
        timetmp = datetime(ts(it),'convertfrom','posixtime');
        pydate = py.datetime.datetime(int32(timetmp.Year),int32(timetmp.Month),int32(timetmp.Day),int32(timetmp.Hour),int32(timetmp.Minute));
        glat = loc(1);
        glon = loc(2);
        fc = py.flipchem.Flipchem(pydate);
        for ih = 1:length(h)
            galt = h(ih);
	    ne = max(pp(ih,it),1e7);
	    ti = max(par(ih,2,it),50);
	    te = max(par(ih,3,it),50);
            outputs = fc.get_point(glat,glon,galt,ne,te,ti);
            outputsm = cell(outputs);
            modelout(ih,9,it) = outputsm{5}; % O2+
            modelout(ih,8,it) = outputsm{6}; % NO+
        end
    end
    
end

