function ElSpecGlowOut = mergeElspecGlowOutputFiles(varargin)
%
% Merge several ElSpec output files.
%
% ElSpecGlowOut = mergeElspecGlowOutputFiles(dfiles)
%
% INPUT:
%  dfiles a comma-separated list of ElSpec-GLOW output files
%
%
% OUTPUT:
%  ElSpecGlowOut a standad ElSpec output list, which is also written
%            to the output file. See ElSpec for details.
%
%
% Details:
%  The merged data are written in a file with the file name
%    ElSpec_<starttime>-<endtime>_merged_<mergetime>.mat
%  For example
%    ElSpec_20180101T012345-20180101T123456_merged_20190204T151617.mat
%
% IV 2019, 2025

% the output struct
ElSpecGlowOut = struct();

% read all the data
ndf = length(varargin);
outlist = cell(ndf);
if ndf==0
return
end
for k=1:ndf
    tmplist = load(varargin{k});
    outlist{k} = tmplist.ElSpecGlowOut;
end
clear tmplist

% try to find common grids for height and energy, and find the
% correct order for the data files
htmp = [];
GLOWhtmp = [];
etmp = [];
tstarts = zeros(ndf,1);
nt = 0;
for k=1:ndf
    htmp = [htmp,outlist{k}.h];
    GLOWhtmp = [GLOWhtmp,outlist{k}.GLOWh];
    etmp = [etmp,outlist{k}.Ec];
    tstarts(k) = outlist{k}.ts(1);
    nt = nt + length(outlist{k}.ts);
end
ElSpecGlowOut.h = unique(htmp);
nh = length(ElSpecGlowOut.h);
ElSpecGlowOut.GLOWh = unique(GLOWhtmp);
nhGLOW = length(ElSpecGlowOut.GLOWh)
ElSpecGlowOut.Ec = unique(etmp);
nE = length(ElSpecGlowOut.Ec);
[dummy,tinds] = sort(tstarts);

% create empty arrays
ElSpecGlowOut.ts = NaN(nt,1);
ElSpecGlowOut.te = NaN(nt,1);
ElSpecGlowOut.pp = NaN(nh,nt);
ElSpecGlowOut.ppstd = NaN(nh,nt);
ElSpecGlowOut.ne = NaN(nh,nt);
ElSpecGlowOut.Ie = NaN(nE,nt);
ElSpecGlowOut.IeStd = NaN(nE,nt);
ElSpecGlowOut.chisqr = NaN(1,nt);
ElSpecGlowOut.FAC = NaN(1,nt);
ElSpecGlowOut.FACstd = NaN(1,nt);
ElSpecGlowOut.E0 = NaN(1,nt);
ElSpecGlowOut.Pe = NaN(1,nt);
ElSpecGlowOut.PeStd = NaN(1,nt);
ElSpecGlowOut.emin = Inf;
ElSpecGlowOut.q = NaN(nh,nt);


ElSpecGlowOut.GLOWtime = NaT(nt,1);
ElSpecGlowOut.GLOWe4278 = NaN(nhGLOW,nt);
ElSpecGlowOut.GLOWe4278max = NaN(1,nt);
ElSpecGlowOut.GLOWh4278max = NaN(1,nt);
ElSpecGlowOut.GLOWh4278mean = NaN(1,nt);
ElSpecGlowOut.GLOWe5577 = NaN(nhGLOW,nt);
ElSpecGlowOut.GLOWe5577max = NaN(1,nt);
ElSpecGlowOut.GLOWh5577max = NaN(1,nt);
ElSpecGlowOut.GLOWh5577mean = NaN(1,nt);
ElSpecGlowOut.GLOWne = NaN(nhGLOW,nt);
ElSpecGlowOut.GLOWnemax = NaN(1,nt);
ElSpecGlowOut.GLOWhnemax = NaN(1,nt);


% interpolate to the common height and energy grids. Throw an error
% if the analysis periods overlap
tcur = 1;
for k=1:ndf

    ii = tinds(k);

    % check that there is no overlap
    if tcur>1
        if outlist{ii}.ts < ElSpecGlowOut.ts(tcur-1)
            error(['The output filest to be merged overlap in time']);
        end
    end

    % add the data to ElSpecGlowOut
    ntcur = length(outlist{ii}.ts);
    tend = tcur+ntcur-1;
    ElSpecGlowOut.ts(tcur:tend) = outlist{ii}.ts;
    ElSpecGlowOut.te(tcur:tend) = outlist{ii}.te;
    ElSpecGlowOut.pp(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.pp,ElSpecGlowOut.h);
    ElSpecGlowOut.ppstd(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.ppstd,ElSpecGlowOut.h);
    ElSpecGlowOut.ne(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.ne,ElSpecGlowOut.h);
    ElSpecGlowOut.Ie(:,tcur:tend) = interp1(outlist{ii}.Ec,outlist{ii}.Ie,ElSpecGlowOut.Ec);
    ElSpecGlowOut.IeStd(:,tcur:tend) = interp1(outlist{ii}.Ec,outlist{ii}.IeStd,ElSpecGlowOut.Ec);
    ElSpecGlowOut.chisqr(tcur:tend) = outlist{ii}.chisqr;
    ElSpecGlowOut.FAC(tcur:tend) = outlist{ii}.FAC;
    ElSpecGlowOut.FACstd(tcur:tend) = outlist{ii}.FACstd;
    ElSpecGlowOut.E0(tcur:tend) = outlist{ii}.E0;
    ElSpecGlowOut.Pe(tcur:tend) = outlist{ii}.Pe;
    ElSpecGlowOut.PeStd(tcur:tend) = outlist{ii}.PeStd;
    ElSpecGlowOut.q(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.q,ElSpecGlowOut.h);

    ElSpecGlowOut.emin = min( ElSpecGlowOut.emin , outlist{ii}.emin );


    ElSpecGlowOut.GLOWtime(tcur:tend) = outlist{ii}.GLOWtime;
    ElSpecGlowOut.GLOWe4278(:,tcur:tend) = outlist{ii}.GLOWe4278;
    ElSpecGlowOut.GLOWe4278max(tcur:tend) = outlist{ii}.GLOWe4278max;
    ElSpecGlowOut.GLOWh4278max(tcur:tend) = outlist{ii}.GLOWh4278max;
    ElSpecGlowOut.GLOWh4278mean(tcur:tend) = outlist{ii}.GLOWh4278mean;
    ElSpecGlowOut.GLOWe5577(:,tcur:tend) = outlist{ii}.GLOWe5577;
    ElSpecGlowOut.GLOWe5577max(tcur:tend) = outlist{ii}.GLOWe5577max;
    ElSpecGlowOut.GLOWh5577max(tcur:tend) = outlist{ii}.GLOWh5577max;
    ElSpecGlowOut.GLOWh5577mean(tcur:tend) = outlist{ii}.GLOWh5577mean;
    ElSpecGlowOut.GLOWne(:,tcur:tend) = outlist{ii}.GLOWne;
    ElSpecGlowOut.GLOWnemax(:,tcur:tend) = outlist{ii}.GLOWnemax;
    ElSpecGlowOut.GLOWhnemax(:,tcur:tend) = outlist{ii}.GLOWhnemax;


    
    tcur = tend + 1;

end

ElSpecGlowOut.radar = outlist{1}.radar;

outfilename = ['ElSpec_',datestr(datetime(round(ElSpecGlowOut.ts(1)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(ElSpecGlowOut.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_merged_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];
save(outfilename,'ElSpecGlowOut','-v7.3');


return
