function ElSpecOut = cutElSpecOutputFile(dfile,starttime,endtime)
%
% Cut a part from an ElSpec output file.
%
% ElSpecOut = cutElSpecOutputFile(dfile,starttime,endtime)
%
% INPUT:
%  dfile      an ElSpec output files
%  starttime  start time of the cut file (datetime)
%  endtime    end time of the cut file (datetime)
%
%
% OUTPUT:
%  ElSpecOut a standad ElSpec output list, which is also written
%            to the output file. See ElSpec for details.
%
%
% Details:
%  The cut data are written in a file with the file name
%    ElSpec_<starttime>-<endtime>_cut_<mergetime>.mat
%  For example
%    ElSpec_20180101T012345-20180101T123456_cut_20190204T151617.mat
%
% IV 2019

% the output struct
ElSpecOut = struct();

% read all the data
tmplist = load(dfile);

ts = posixtime(starttime);
te = posixtime(endtime);
iit = tmplist.ElSpecOut.ts >= ts & tmplist.ElSpecOut.te <= te;


% create empty arrays
ElSpecOut.h = tmplist.ElSpecOut.h;
ElSpecOut.Ec = tmplist.ElSpecOut.Ec;
ElSpecOut.dE = tmplist.ElSpecOut.dE;
ElSpecOut.ts = tmplist.ElSpecOut.ts(iit);
ElSpecOut.te = tmplist.ElSpecOut.te(iit);
ElSpecOut.pp = tmplist.ElSpecOut.pp(:,iit);
ElSpecOut.ppstd = tmplist.ElSpecOut.ppstd(:,iit);
ElSpecOut.ne = tmplist.ElSpecOut.ne(:,iit);
ElSpecOut.Ie = tmplist.ElSpecOut.Ie(:,iit);
ElSpecOut.IeStd = tmplist.ElSpecOut.IeStd(:,iit);
ElSpecOut.chisqr = tmplist.ElSpecOut.chisqr(iit);
ElSpecOut.FAC = tmplist.ElSpecOut.FAC(iit);
ElSpecOut.FACstd = tmplist.ElSpecOut.FACstd(iit);
ElSpecOut.E0 = tmplist.ElSpecOut.E0(iit);
ElSpecOut.Pe = tmplist.ElSpecOut.Pe(iit);
ElSpecOut.PeStd = tmplist.ElSpecOut.PeStd(iit);
ElSpecOut.emin = tmplist.ElSpecOut.emin;;
ElSpecOut.q = tmplist.ElSpecOut.q(:,iit);;

ElSpecOut.radar = tmplist.ElSpecOut.radar;

outfilename = ['ElSpec_',datestr(datetime(round(ElSpecOut.ts(1)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(ElSpecOut.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_cut_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];
save(outfilename,'ElSpecOut','-v7.3');


return
