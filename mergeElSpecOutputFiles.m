function ElSpecOut = mergeElSpecOutputFiles(varargin)
%
% Merge several ElSpec output files.
%
% ElSpecOut = mergeOutputFiles(dfiles)
%
% INPUT:
%  dfiles a comma-separated list of ElSpec output files
%
%
% OUTPUT:
%  ElSpecOut a standad ElSpec output list, which is also written
%            to the output file. See ElSpec for details.
%
%
% Details:
%  The merged data are written in a file with the file name
%    ElSpec_<starttime>-<endtime>_merged_<mergetime>.mat
%  For example
%    ElSpec_20180101T012345-20180101T123456_merged_20190204T151617.mat
%
% IV 2019

% the output struct
ElSpecOut = struct();

% read all the data
ndf = length(varargin);
outlist = cell(ndf);
if ndf==0
return
end
for k=1:ndf
    tmplist = load(varargin{k});
    outlist{k} = tmplist.ElSpecOut;
end
clear tmplist

% try to find common grids for height and energy, and find the
% correct order for the data files
htmp = [];
etmp = [];
tstarts = zeros(ndf,1);
nt = 0;
for k=1:ndf
    htmp = [htmp,outlist{k}.h];
    etmp = [etmp,outlist{k}.Ec];
    tstarts(k) = outlist{k}.ts(1);
    nt = nt + length(outlist{k}.ts);
end
ElSpecOut.h = unique(htmp);
nh = length(ElSpecOut.h);
ElSpecOut.Ec = unique(etmp);
nE = length(ElSpecOut.Ec);
[dummy,tinds] = sort(tstarts);

% create empty arrays
ElSpecOut.ts = NaN(nt,1);
ElSpecOut.te = NaN(nt,1);
ElSpecOut.pp = NaN(nh,nt);
ElSpecOut.ppstd = NaN(nh,nt);
ElSpecOut.ne = NaN(nh,nt);
ElSpecOut.Ie = NaN(nE,nt);
ElSpecOut.IeStd = NaN(nE,nt);
ElSpecOut.chisqr = NaN(1,nt);
ElSpecOut.FAC = NaN(1,nt);
ElSpecOut.FACstd = NaN(1,nt);
ElSpecOut.Pe = NaN(1,nt);
ElSpecOut.PeStd = NaN(1,nt);
ElSpecOut.emin = Inf;
ElSpecOut.q = NaN(nh,nt);

% interpolate to the common height and energy grids. Throw an error
% if the analysis periods overlap
tcur = 1;
for k=1:ndf

    ii = tinds(k);

    % check that there is no overlap
    if tcur>1
        if outlist{ii}.ts < ElSpecOut.ts(tcur-1)
            error(['The output filest to be merged overlap in ' ...
                   'time']);
        end
    end

    % add the data to ElSpecOut
    ntcur = length(outlist{ii}.ts);
    tend = tcur+ntcur-1;
    ElSpecOut.ts(tcur:tend) = outlist{ii}.ts;
    ElSpecOut.te(tcur:tend) = outlist{ii}.te;
    ElSpecOut.pp(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.pp,ElSpecOut.h);
    ElSpecOut.ppstd(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.ppstd,ElSpecOut.h);
    ElSpecOut.ne(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.ne,ElSpecOut.h);
    ElSpecOut.Ie(:,tcur:tend) = interp1(outlist{ii}.Ec,outlist{ii}.Ie,ElSpecOut.Ec);
    ElSpecOut.IeStd(:,tcur:tend) = interp1(outlist{ii}.Ec,outlist{ii}.IeStd,ElSpecOut.Ec);
    ElSpecOut.chisqr(tcur:tend) = outlist{ii}.chisqr;
    ElSpecOut.FAC(tcur:tend) = outlist{ii}.FAC;
    ElSpecOut.FACstd(tcur:tend) = outlist{ii}.FACstd;
    ElSpecOut.Pe(tcur:tend) = outlist{ii}.Pe;
    ElSpecOut.PeStd(tcur:tend) = outlist{ii}.PeStd;
    ElSpecOut.q(:,tcur:tend) = interp1(outlist{ii}.h,outlist{ii}.q,ElSpecOut.h);

    ElSpecOut.emin = min( ElSpecOut.emin , outlist{ii}.emin );

    tcur = tend + 1;

end

outfilename = ['ElSpec_',datestr(datetime(round(ElSpecOut.ts(1)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(ElSpecOut.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_merged_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];
save(outfilename,'ElSpecOut','-v7.3');


return