function ElSpecOut = ElSpec(varargin)
%
% ElSpec Estimate primary energy spectra of precipitating electrons from
% EISCAT incoherent scatter data.
%
%
% The primary energy spectra are modeled  using a parametrix model
% of the form
% s = E.*exp(P_n(E))
% where P_n is an n'th order polynomial.
%
% Corrected Akaike information criterion is used for seleting the
% optimal order of the polynomial P_n.
%
% ElSpecOut = ElSpec(varargin)
%
%
% INPUT:
%
%  name-value pairs, the available arguments are listed below.
%
%  At least  ppdir OR fitdir must be defined.
%
%  For example ElSpec('ppdir',<path-to-pp-files>)
%
%  ppdir        Path to GUISDAP raw densities. See details.
%  fitdir       Path to GUISDAP parameter fit results. See details.
%  experiment   EISCAT experiment name, e.g. 'beata' or 'arc1'
%  radar        Radar, 'uhf', 'vhf', or 'esr'
%  version      Experiment version number.
%  hmin         Minimum altitude in km, must be between 80 and 150 km.
%  hmax         Maximum altitude in km, must be between 80 and 150 km.
%  btime        Analysis start time [year,month,day,hour,minute,second]
%  etime        Analysis end time [year,month,day,hour,minute,second]
%  ionomodel    Ion production model, 'Fang' or 'Sergienko'. See details.
%  recombmodel  Recombination model, 'Rees', 'ReesO2+', 'ReesN2+', 'ReesNO+',
%               'SheehanGrO2+', 'SheehanGrN2+', 'SheehanGrNO+', 'SheehanGrFlipchem'
%               'delPozo1', or 'delPozo2'. See details.
%  integtype   Type of continuity function integration, 'endNe',
%              'integrate', or 'equilibrium'. See details.
%  egrid        Energy grid [eV]
%  maxorder     range of the number of nodes in the spline
%               expansion.
%  plotres      plot results with ElSpecPlot? 1=plot, 0=no plot,
%               default 1
%  tres         "type" of time resolution 'best' or 'dump' [default]
%  emin         smallest energy to include in FAC and power
%               estimation, default 1 keV
%  ieprior      prior model test....,requires a lot of manual
%                work... use empty array to skip the prior...
%  stdprior
%  ninteg       number of ne-profiles to use in each integration
%               period, default 6
%  nstep        number of ne-slices in each time-step, default 1
%  saveiecov    logical, should the large covariance matrices of
%               the flux estimates be saved? default false.
%  fadev        maximum beam direction deviation from field-aligned  [deg], default 3
%  bottomstdfail standard deviation given for model Ne in the lowest observed gate when the
%                guisdap fit has failed [m^-3], default 1e10
%  refineEgrid  logial, refine energy grid according to the lowest measured altitude? Default 1
%
% OUTPUT:
%  ElSpecOut    A MATLAB structure with fields:...
%
%
%
% Details:
%
%  The analysis uses two sets of GUISDAP results as input: a set of
%  high-time resolution raw electron densities, and a set of lower
%  resolution 'normal' fit results. These cannot be produced in a
%  single GUISDAP run, because time resolution of the raw densities
%  would then be matched to time resolution of the full fit. The
%  raw densities in ppdir can be produced by means of setting
%
%  analysis_ppshortlags=1
%
%  and
%
%  analysis_pponly=1
%
%  in GUISDAP. The
%  GUISDAP output will be experiment- and
%  radar-specific. Currently, we have a readers for uhf beata and
%  arc1 only.
%
%  Two models for ion production profiles by monoenergetic
%  electrons are implemented. These are Sergienko1993 and
%  Fang2010. See
%  help ion_production_Sergienko1993
%  and
%  help ion_production_Fang2010
%  for details.
%
%  Nine recombination models are implemented. 'Rees' is based on
%  ion abundances from IRI and values from Rees
%  (1989). 'SheehanGr' uses values for ground state ions from
%  Sheehan and St.-Maurice (2004). 'SheehanEx' uses values for
%  vibrationally excited ions from Sheehan and St.-Maurice (2004).
%  'ReesO2+','ReesN2+', and 'ReesNO+' are Rees recombination rates for
%  pure O2+, N2+, adn NO+, correspondingly. 'SheehanGrO2+',
%  'SheehanGrN2+' and 'SheehanGrNO+' are the corresponding Sheehan
%  and St.-Maurice (2004) recombination rates for vibrational
%  ground-states. 'SheehanGrFlipchem' is Sheehan and St.-Maurice (2004)
%  with ion composition from the flipchem model.
%  'delPozo1' and 'delPozo2' are simple analytic approximations. See
%   help effective_recombination_rate
%  for details.
%
%  Three different methods for integrating the electron continuity
%  equation are available. 'endNe' ignores the actual integration,
%  and returns the electron density profile at end of the time
%  step. 'integrate' returns average of the time-dependent density
%  profile over the time-step dt. 'equilibrium' ignores the initial
%  condition and assumes a photochemical equilibrium. 'equilibrium'
%  is always forced for the very first time step, where the initial
%  condition is inknown. See
%  help integrate_continuity_equation
%
%
%  The software uses tabulated values of IRI parameters and uses
%  the MATLAB Aerospace toolbox implementation of MSIS model. The
%  IRI tables are provided as part of the distribution,
%  and are automatically read by the function 'modelParams'. An R
%  function for creating new tables, createModelParameterFile, is
%  provided within the MATLAB distribution. The R function uses
%  packages IRI2016 and R.matlab. The latter is a standard CRAN
%  pacakge, whereas the package IRI2016 can be requested from Ilkka
%  Virtanen (ilkka.i.virtanen@oulu.fi).
%
%
%  The bottomstdfail argument is used to suppress unrealistic high-energy
%  peaks from the energy spectra when GUISDAP fits in the bottom part of the
%  Ne profile have failed. It is usually safe to assume a very small electron density
%  with a small variance, because the most common reason for failed fits in the
%  bottom part is very low Ne. This parameter has no effect on analysis with
%  power profile input, because there are no failed fits GUISDAP raw density data. 
%
% IV 2017, 2018, 2021
%
% Copyright I Virtanen <ilkka.i.virtanen@oulu.fi> and B Gustavsson <bjorn.gustavsson@uit.no>
% This is free software, licensed under GNU GPL version 2 or later

persistent ElSpecLicenseShown




% display the copyright message only once
if isempty(ElSpecLicenseShown)
    ElSpecLicenseShown = true;
    disp('  ')
    disp('Copyright 2015-2018, Ilkka Virtanen and Bjorn Gustavsson')
    disp(['This is free software, licensed under GNU GPL version 2 or later.'])
    disp('  ')
    disp(['This program is distributed in the hope that it will be ' ...
          'useful, '])
    disp(['but WITHOUT ANY WARRANTY; without even the implied ' ...
          'warranty of '])
    disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.')
    disp('See the GNU General Public License for more details')
    disp('  ')
    disp('  ')
    disp('  ')
    disp('Reference: I.I. Virtanen, B. Gustavsson, et al., Electron ')
    disp('energy spectrum and auroral power estimation from ')
    disp('incoherent scatter radar measurements, J. Geophys. Res. ')
    disp('Space Physics, 123, 6865–6887, 2018, doi:10.1029/2018JA025636')
    disp('  ')
    pause(5)
end

% parse the input
p = inputParser;

% power profiles
defaultPPdir = [];
checkPPdir = @(x) (exist(x,'dir')|isempty(x)|strcmp(x,'simu'));

% plasma parameters
defaultFitdir = [];
checkFitdir = @(x) (exist(x,'dir')|isempty(x)|(exist(x,'file')&strcmp( x((end-4):end) ,'.hdf5'))|(exist(x,'file')&strcmp( x((end-3):end) ,'.mat')));

% name of the experiment
defaultExperiment = 'beata';
checkExperiment = @(x) (isstr(x));

% radar, only uhf exps implemented at the moment
defaultRadar = 'uhf';
validRadar = {'uhf','vhf','esr','42m','32m','esr42m','esr32m'};
checkRadar = @(x) any(validatestring(x,validRadar));

% experiment version
defaultVersion = 1;
checkVersion = @(x) (isnumeric(x) & length(x)==1);

% minimum altitude
defaultHmin = 80;
checkHmin = @(x) (x>=80 & x <150);

% maximum altitude
defaultHmax = 150;
checkHmax = @(x) ( x>80 & x<=150);

% ion production model
defaultIonomodel = 'Fang';
validIonomodel = {'Fang','Sergienko','Rees'};
checkIonomodel = @(x) any(validatestring(x,validIonomodel));

% recombination model
defaultRecombmodel = 'SheehanGr';
validRecombmodel = {'SheehanGr','SheehanEx','Rees','ReesO2+', ...
                    'ReesN2+','ReesNO+','SheehanGrO2+', ...
                    'SheehanGrN2+','SheehanGrNO+','SheehanGrFlipchem',...
                    'delPozo1','delPozo2'};
checkRecombmodel = @(x) any(validatestring(x,validRecombmodel));

% type of integration
defaultIntegtype = 'integrate';
validIntegtype = {'integrate','equilibrium','endNe'};
checkIntegtype = @(x) any(validatestring(x,validIntegtype));

% energy grid
defaultE = logspace(1,5,300);
checkE = @(x) (isnumeric(x) & length(x)>1 & all(x>0));

% maximum order of the polynomial model
defaultMaxorder = 5;
checkMaxorder = @(x) (length(x)==1 & all(x>=1));

% plot the results?
defaultPlotres = true;
checkPlotres = @(x) ((isnumeric(x)|islogical(x))&length(x)==1);

% time resolution (best available or EISCAT dump length)
defaultTres = 'dump';
validTres = {'dump','best'};
checkTres = @(x) any(validatestring(x,validTres));

% minimum energy for FAC and power calculation
defaultEmin = 1e3;
checkEmin = @(x) (isnumeric(x) & length(x)==1 & x>0);

% a priori differential number flux
defaultIeprior = [];
checkIeprior = @(x) (ismatrix(x));

% a priori standard deviations
defaultStdprior = [];
checkStdprior = @(x) (ismatrix(x));

% number of time-slices integrated together
defaultNinteg = 6; % 30 s in beata, note that only the first is used
                   % with full weight. the time resolution will be
                   % effectively one slice for high Ne
checkNinteg = @(x) (isnumeric(x) & length(x)==1 & all(x>0));

% number of time-slices to proceed on each time-step
defaultNstep = 1;
checkNstep = @(x) (isnumeric(x) & length(x)==1 & all(x>0));

% save the  large flux covariance matrix?
defaultSaveiecov = false;
checkSaveiecov = @(x) ((isnumeric(x)|islogial(x))&length(x)==1);

% use only field-aligned data?
defaultFAdev = 3;
checkFAdev = @(x) (isnumeric(x)&length(x)==1);

% standard deviation for model Ne in the lowest observed altitude when the guisdap fit has failed
defaultBottomStdFail = 1e10;
checkBottomStdFail = @(x) (isnumeric(x)&length(x)==1 & x>0);

% should the energy axis be refined according to the lowest altitude=
defaultRefineEgrid = true;
checkRefineEgrid = @(x) ((isnumeric(x)|islogial(x))&length(x)==1);

% start time
defaultBtime = [];

% end time
defaultEtime = [];

% function for checking the time limits
checkTlim = @(x) (isnumeric(x) & length(x)==6 & all(x>=0));

addParameter(p,'btime',defaultBtime,checkTlim);
addParameter(p,'etime',defaultEtime,checkTlim);
addParameter(p,'ppdir',defaultPPdir,checkPPdir);
addParameter(p,'fitdir',defaultFitdir,checkFitdir);
addParameter(p,'experiment',defaultExperiment,checkExperiment);
addParameter(p,'radar',defaultRadar,checkRadar);
addParameter(p,'version',defaultVersion,checkVersion);
addParameter(p,'hmin',defaultHmin,checkHmin);
addParameter(p,'hmax',defaultHmax,checkHmax);
addParameter(p,'ionomodel',defaultIonomodel,checkIonomodel);
addParameter(p,'recombmodel',defaultRecombmodel,checkRecombmodel);
addParameter(p,'integtype',defaultIntegtype,checkIntegtype);
addParameter(p,'egrid',defaultE,checkE);
addParameter(p,'maxorder',defaultMaxorder,checkMaxorder);
addParameter(p,'plotres',defaultPlotres,checkPlotres);
addParameter(p,'tres',defaultTres,checkTres);
addParameter(p,'emin',defaultEmin,checkEmin);
addParameter(p,'ieprior',defaultIeprior,checkIeprior);
addParameter(p,'stdprior',defaultStdprior,checkStdprior);
addParameter(p,'ninteg',defaultNinteg,checkNinteg);
addParameter(p,'nstep',defaultNstep,checkNstep);
addParameter(p,'saveiecov',defaultSaveiecov,checkSaveiecov);
addParameter(p,'fadev',defaultFAdev,checkFAdev);
addParameter(p,'bottomstdfail',defaultBottomStdFail,checkBottomStdFail);
addParameter(p,'refineEgrid',defaultRefineEgrid,checkRefineEgrid);
parse(p,varargin{:})

%out = struct;

% collect all the inputs to the output list
out = p.Results;
out.E = out.egrid;

% check that we have either ppdir or fitdir, or both
if isempty(out.ppdir) & isempty(out.fitdir)
    error('Either ppdir or fitdir must be given')
end

% we have input routines only for uhf arc1 and beata power profiles
if ~isempty(out.ppdir)
    wrongexp = false;
    if strcmp(out.radar,'uhf')
        if ~any(strcmp(out.experiment,{'beata','arc1'}))
            wrongexp = true;
        end
    else
        wrongexp = true;
    end
    if wrongexp
        error(['Power profile input routines are available only for ' ...
               'UHF arc1 and UHF beata experiments.'])
    end
end

% check that nstep is not larger than ninteg
if out.nstep > out.ninteg
    error('nstep must not be larger than ninteg.')
end



% read the data and model values
% a hack to allow reading in simulated data
if strcmp(out.ppdir,'simu')
    load(out.fitdir);
    out.h = simudata.h;
    out.ts = simudata.ts;
    out.te = simudata.te;
    out.pp = simudata.pp;
    out.ppstd = simudata.ppstd;
    out.par = simudata.par;
    out.parstd = simudata.parstd;
    out.iri = simudata.iri;
else
    readIRI = false;
    if any(strcmp(out.recombmodel,{'SheehanGr','SheehanEx','Rees'}))
        readIRI = true;
    end
    [out.h,out.ts,out.te,out.pp,out.ppstd,out.par,out.parstd,out.iri,out.f107,out.f107a,out.f107p,out.ap,out.loc,out.azel,out.I] = ...
        readFitData( out.ppdir , out.fitdir , out.hmin , out.hmax , ...
                     out.btime , out.etime , out.experiment , out.radar , ...
                     out.version , out.tres , readIRI, p.Results.fadev , p.Results.bottomstdfail);
    if strcmp(out.recombmodel,'SheehanGrFlipchem')
    out.iri = calculateFlipchemComposition(out.ts,out.h,out.par,out.pp,out.loc,out.iri);
    end
    % % warn about the ESR compositions
    % if strcmp(p.Results.radar,'esr')
    %     if readIRI
    %         disp(['WARNING: With the selected recombination model ElSpec calculates the ion  compositions ' ...
    %               'with EISCAT Tromso coordinates, the compositions may be ' ...
    %               'incorrect for the ESR radar'])
    %         % disp(['WARNING: ElSpec calculates the ion and neutral compositions ' ...
    %         %       'with EISCAT Tromso coordinates, the compositions may be ' ...
    %         %       'incorrect for the ESR radar'])
    %     end
    % end


end
if isempty(out.h)
    disp('No data')
    return
end
% disp('************************************')
% disp('TESTING with SIC composition, ElSpec, line 318!')
% disp('************************************')
% %load('/Users/ilkkavir/Documents/EISCAT_LAB/ElectronPrecipitationSpectra/SIC_composition/Ilkalle_composition.mat');
% load('Ilkalle_composition.mat');
% tsic = posixtime(datetime(t,'convertfrom','datenum'));
% %interpolate to the correct height grid
% O2p = zeros(length(out.h),length(tsic));
% Op = zeros(length(out.h),length(tsic));
% NOp = zeros(length(out.h),length(tsic));
% for indd = 1:length(tsic)
%    O2p(:,indd) = interp1(Hgt,O2plus(:,indd),out.h,'linear','extrap');
%    Op(:,indd) = interp1(Hgt,Oplus(:,indd),out.h,'linear','extrap');
%    NOp(:,indd) = interp1(Hgt,NOplus(:,indd),out.h,'linear','extrap');
% end
% %interplation in time
% for indd = 1:length(out.h)
%    out.iri(indd,8,:) = interp1(tsic,NOp(indd,:),out.ts,'linear','extrap');
%    out.iri(indd,9,:) = interp1(tsic,O2p(indd,:),out.ts,'linear','extrap');
%    out.iri(indd,10,:) = interp1(tsic,Op(indd,:),out.ts,'linear','extrap');
% end



if out.refineEgrid
    
    % Calculate ion production matrix to approximately adjust the energy grid
    [A,Ec,dE] = ion_production( out.E , out.h*1000 , out.iri(:,4,1) , ...
                                out.iri(:,5,1) , out.iri(:,6,1) , out.iri(:,7,1) , ...
                                out.iri(:,1,1) , out.ionomodel , out.I);
    % adjust the energy grid according to the lowest measured altitude
    [~,imax] = max(A(1,:));
    if imax < length(Ec)
        Ec = Ec(1:imax);
        dE = dE(1:imax);
        A = A(:,1:imax);
        out.E = out.E(1:(imax+1));
        out.egrid = out.egrid(1:(imax+1));
        %    disp([imax length(Ec)])
    end

end

% time step sizes
if length(out.te)==1
    dt = out.te-out.ts(1);
else
    out.dt = diff( out.te );
    out.dt = [out.dt(1);out.dt(:)];
end

% some dimensions and initializations
nt = length(out.ts); % number of time steps
nh = length(out.h);  % number of heights
nE = length(out.E) - 1; % number of energy bins
out.ne = NaN(nh,nt);   % an array for electron density estimates
out.neEnd = NaN(nh,nt);   % an array for electron density estimates
                          % at integration end points
out.neEndStd = NaN(nh,nt);   % an array for electron density estimates
out.neEndCov = NaN(nh,nh,nt);   % an array for electron density estimates

out.Ie = NaN(nE,nt);   % an array for flux estimates
out.AICc = NaN(out.maxorder,nt); % an array for the information
                                  % criterion values

out.IeCov = NaN(nE,nE,nt);   % Covariances of flux estimates
out.IeStd = NaN(nE,nt); % standard deviation of flux estimates
out.alpha = NaN(nh,nt); % an array for the effective recombination
                        % rates
out.q = NaN(nh,nt); % an array for ion production rates
out.polycoefs = NaN(out.maxorder,out.maxorder+1,nt);
out.polycoefsCovar = NaN(out.maxorder,out.maxorder+1,out.maxorder+1,nt);
out.exitflag = NaN(out.maxorder,nt);
out.output = cell(out.maxorder,nt);
out.best_order = NaN(1,nt);
out.chisqr = NaN(1,nt);
out.FAC = NaN(1,nt);
out.FACstd = NaN(1,nt);
out.Pe = NaN(1,nt);
out.PeStd = NaN(1,nt);

% differential number flux, differential energy flux, and peak energy
% Nlux will be a copy of Ie, but also Ie is kept for backward compatibility
out.Nflux = NaN(nE,nt);
out.Eflux = NaN(nE,nt);
out.E0 = NaN(1,nt);

% options for the MATLAB fit routines
fms_opts = optimset('fminsearch');
fms_opts.Display = 'off';%'off';%'final';
fms_opts.MaxFunEvals=1e4;
fms_opts.MaxIter=1e6;
fms_opts.TolFun=1e-8;
fms_opts.TolX=1e-8;

% name of the output file
out.outputfile = ['ElSpec_',datestr(datetime(round(out.ts(1)), ...
                                            'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'-',datestr(datetime(round(out.te(end)),'ConvertFrom','posixtime'),'yyyymmddTHHMMss'),'_',out.experiment,'_',out.radar,'_',out.ionomodel,'_',out.recombmodel,'_',out.integtype,'_',num2str(out.ninteg),'_',num2str(out.nstep),'_',out.tres,'_',datestr(datetime('now'),'yyyymmddTHHMMSS'),'.mat'];

% Then print all the inputs
fprintf('\nElSpec input arguments\n\n')

fprintf('%12s: %s\n','ppdir',out.ppdir);
fprintf('%12s: %s\n','fitdir',out.fitdir);
fprintf('%12s: %s\n','experiment',out.experiment);
fprintf('%12s: %s\n','radar',out.radar);
fprintf('%12s: %i\n','version',out.version);
fprintf('%12s: %5.1f\n','hmin [km]',out.hmin);
fprintf('%12s: %5.1f\n','hmax [km]',out.hmax);
fprintf('%12s: ','btime');
for kk=1:length(out.btime)
    fprintf('%i ',out.btime(kk));
end
fprintf('\n');
fprintf('%12s: ','etime');
for kk=1:length(out.etime)
    fprintf('%i ',out.etime(kk));
end
fprintf('\n');
fprintf('%12s: %s\n','ionomodel',out.ionomodel);
fprintf('%12s: %s\n','recombmodel',out.recombmodel);
fprintf('%12s: %s\n','integtype',out.integtype);
fprintf('%12s: ','E [keV]');
for kk=1:(length(out.egrid)/10)
    if kk>1
        fprintf('%14s','');
    end
    fprintf(' %5.2f',out.E((kk-1)*10+1:min(length(out.E),kk*10))/1000);
    fprintf('\n');
end
fprintf('%12s: %i\n','maxorder',out.maxorder);
fprintf('%12s: %i\n','plotres',out.plotres);
fprintf('%12s: %s\n','tres',out.tres);
fprintf('%12s: %10.0f\n','Emin [eV]',out.emin);
fprintf('%12s: %i\n','ninteg',out.ninteg);
fprintf('%12s: %i\n','nstep',out.nstep);

fprintf('Will write results in:\n %s\n',out.outputfile);

% save interval is 100 step, independently from the time resolution
ndtsave = 100;%ceil(mean(120./out.dt));

% iteration over all time steps
for tt = 1:out.nstep:nt-out.ninteg

    % initialize temporary variables for ne and Ie
    ne_tmp = NaN(nh,out.maxorder);
    IE_tmp = NaN(nE,out.maxorder);

    % Update the ion production matrix
    [A,Ec,dE] = ion_production( out.E , out.h*1000 , out.iri(:,4,tt) , ...
                                out.iri(:,5,tt) , out.iri(:,6,tt) , out.iri(:,7,tt) , ...
                                out.iri(:,1,tt) , out.ionomodel , out.I);
    A(isnan(A)) = 0;

    % save an example of A on out
    if tt==1
        out.A = A;
    end

    % update the effective recombination coefficient. Assume N2+
    % density to be zero
    out.alpha(:,tt) = effective_recombination_coefficient( out.h , ...
                         out.par(:,3,tt) , ...
                         out.iri(:,9,tt) , ...
                         out.iri(:,9,tt).*0 , ...
                         out.iri(:,8,tt) , ...
                         out.recombmodel );
    for istep = 0:out.nstep-1
        out.alpha(:,tt+istep) = out.alpha(:,tt);
    end


    % Initial values, different for the first step and the others
    if tt == 1
        % Use 'equilibrium' as 'integ_type' on the first time step
        integ_type_tt = 'equilibrium';
        ne0 = zeros(nh,1);
        ne0Std = zeros(nh,1)+1e12;
        ne0Cov = diag(ones(nh,1)*1e24);
        % copy Ec and dE to out
        out.Ec = Ec;
        out.dE = dE;

    else % tt > 1
        integ_type_tt = out.integtype;
        ne0 = out.neEnd(:,tt-1);
        ne0Std = out.neEndStd(:,tt-1);
        ne0Cov = out.neEndCov(:,:,tt-1);
    end
    % if the previous error estimation failed, replace ne0Std with
    % large values
    ne0diag = false;
    if any(isnan(ne0Std)) | any(imag(ne0Std)~=0) | any(diag(ne0Cov)<=0)
        ne0Std = zeros(nh,1) + 1e12;
        ne0Cov = diag(ones(nh,1)*1e24);
        ne0diag = true;
        %        disp('replacing ne0Cov with a diagonal matrix')
    end

    % recombination timescale
    tscale = 1./(out.alpha(:,tt).*(ne0(:)+1));
    tscale(tscale>300) = 300;

    % normalize variances
    varscale = out.ppstd(:,tt:tt+out.ninteg-1).*0;
    ttimes = cumsum(out.dt(tt:tt+out.ninteg-1)) - out.dt(tt);
    for indtt=1:out.ninteg
        varscale(:,indtt) = exp( ttimes(indtt) ./ tscale * 4) ;
    end
    if all(ne0==0)
        varscale = varscale*0+1e3;
        varscale(:,1) = 1;
    end
    stdppfit = out.ppstd(:,tt:tt+out.ninteg-1).*sqrt(varscale);

    % effective number of measurements
    nmeaseff = sum(sum(1./varscale));


    % then run the fit for all different numbers of nodes
    for nn = 1:out.maxorder

        % an initial guess of the coefficients
        if nn==1
            X0 = [20,-5];
        else
            X0 = [out.polycoefs(nn-1,1:nn,tt),0];
        end
        % do not limit the coefficients, this is handled with a
        % prior now
        %X0min = -Inf*X0;
        %X0max = Inf*X0;

        % an optioal prior distribution from a steady-state fit...
        if ~isempty(out.ieprior) & ~isempty(out.stdprior)
            [x,fval,exitflag,output] = fminsearch( @(p) ElSpec_fitfun( ...
                p , out.pp(:,tt:tt+out.ninteg-1) , stdppfit , ne0 , A , ...
                out.alpha(:,tt) , out.dt(tt:tt+out.ninteg-1) , out.Ec , ...
                out.dE , integ_type_tt , out.ieprior(:,tt) , out.stdprior(:,tt) ...
                , nmeaseff) , X0 , fms_opts);
        % the normal fit
        else
            [x,fval,exitflag,output] = fminsearch( @(p) ElSpec_fitfun( ...
                p , out.pp(:,tt:tt+out.ninteg-1) , stdppfit , ne0 , A , ...
                out.alpha(:,tt) , out.dt(tt:tt+out.ninteg-1) , out.Ec , ...
                out.dE , integ_type_tt , [] , [] , nmeaseff ...
                ) , X0 , fms_opts);
        end

        % collect results
        for istep = 0:out.nstep-1
            out.polycoefs(nn,1:nn+1,tt+istep) = x;
            out.exitflag(nn,tt+istep) = exitflag;
            out.output{nn,tt+istep} = output;
            out.AICc(nn,tt+istep) = fval;
        end

    end

    % pick the best AICc value
    [minAIC,locminAIC] = min(out.AICc(:,tt));
    out.best_order(tt:tt+out.nstep-1) = locminAIC;


    % error estimation for the coefficients from linearized theory
    if ~isempty(out.ieprior) & ~isempty(out.stdprior)
        xCovar = ElSpec_error_estimate_polycoefs( out.pp(:,tt:tt+out.ninteg-1) ...
                                                  , stdppfit , ne0 , ...
                                                  ne0Cov , A, out.alpha(:,tt) ...
                                                  , out.dt(tt:tt+out.ninteg-1) ...
                                                  , out.Ec , out.dE , ...
                                                  integ_type_tt , out.ieprior(:,tt) ...
                                                  , out.stdprior(:,tt) , out.polycoefs( ...
                                                      locminAIC , ...
                                                      1:locminAIC+1 , ...
                                                      tt ) , nmeaseff ...
                                                  );
    else
        xCovar = ElSpec_error_estimate_polycoefs( out.pp(:,tt:tt+out.ninteg-1) ...
                                                  , stdppfit , ne0 , ...
                                                  ne0Cov , A, out.alpha(:,tt) ...
                                                  , out.dt(tt:tt+out.ninteg-1) ...
                                                  , out.Ec , out.dE , ...
                                                  integ_type_tt , [] ...
                                                  , [] , out.polycoefs( ...
                                                      locminAIC , ...
                                                      1:locminAIC+1 , ...
                                                      tt ) , nmeaseff ...
                                                  );
    end        

    % flux covariance matrix
    out.IeCov(:,:,tt) = ElSpec_error_estimate_Ie( xCovar(1:locminAIC+1 ...
                                                      , 1:locminAIC+1) ...
                                                  , out.polycoefs(locminAIC ...
                                                      , 1:locminAIC+1,tt) ...
                                                  , out.Ec );

    % calculate the flux and electron density profile corresponding
    % the best fit
    for istep = 0:out.nstep-1
        out.polycoefsCovar(locminAIC,1:locminAIC+1,1:locminAIC+1,tt+istep) ...
            = xCovar(1:locminAIC+1,1:locminAIC+1);

        out.IeCov(:,:,tt+istep) = out.IeCov(:,:,tt);
        out.IeStd(:,tt+istep) = sqrt(diag(squeeze(out.IeCov(:,:,tt))));
        out.Ie(:,tt+istep) = model_spectrum( out.polycoefs( locminAIC ...
                                                          , 1:locminAIC+1 ...
                                                          , tt+istep ...
                                                          ) , Ec );
        % correct integration type for the endNe estimation
        switch lower(integ_type_tt)
          case 'endne'
            integ_type_end = 'endNe';
          case 'integrate'
            integ_type_end = 'endNe';
          case 'equilibrium'
            integ_type_end = 'equilibrium';
          case 'linearend'
            integ_type_end = 'linearend';
          case 'linearint'
            integ_type_end = 'linearend';
          otherwise
            error(['Error. \nUnknown integ_type %s in ElSpec' ] , type )
        end

        % modeled Ne at end of the time step
        out.neEnd(:,tt+istep) = integrate_continuity_equation( ne0 , out.Ie(:,tt) ...
                                                          , sum(out.dt(tt:tt+istep)) , A ...
                                                          , out.dE(:) ...
                                                          , out.alpha(:,tt+istep) ...
                                                          , integ_type_end ...
                                                          );

        % standard deviation and covariance matrix fo neEnd
        [out.neEndStd(:,tt+istep) , out.neEndCov(:,:,tt+istep)]= ...
            ElSpec_error_estimate_ne2( ne0 , ...
                                     out.polycoefs(locminAIC , 1:locminAIC+1,tt)  , ...
                                     xCovar, ...
                                     sum(out.dt(tt:tt+istep)) , A , out.Ec , ...
                                     out.dE(:) , out.alpha(:,tt+istep) , ...
                                     integ_type_end );

        % pick the covariance for the starting point of the next integration
        if tt==1
            ne00 = ne0;
            ne00Std = ne0Std;
            ne00Cov = ne0Cov;
        else
            ne00 = out.neEnd(:,tt+istep-1);
            ne00Std = out.neEndStd(:,tt+istep-1);
            ne00Cov = out.neEndCov(:,:,tt+istep-1);
            if ne0diag
                ne00Std = zeros(nh,1)+1e12;
                ne00Cov = diag(ones(nh,1)*1e24);
            end
        end

        % the modeled Ne (which should match with the measurements)
        out.ne(:,tt+istep) = integrate_continuity_equation( ne00 , out.Ie(:,tt) ...
                                                          , out.dt(tt+istep) , A ...
                                                          , out.dE(:) ...
                                                          , out.alpha(:,tt+istep) ...
                                                          , integ_type_tt ...
                                                          );

        % ion production...
        out.q(:,tt+istep) = A*(out.Ie(:,tt+istep).*out.dE');

        %chi-squared
        out.chisqr(tt+istep) = mean(((out.ne(:,tt+istep)-out.pp(:,tt+istep))./out.ppstd(:,tt+istep)).^2);

        % integrate the spectra and multiply with e to get the current
        % density
        Eind_fac = out.Ec >= out.emin;
        out.FAC(tt+istep) = sum(out.Ie(Eind_fac,tt+istep).*out.dE(Eind_fac)')*1.60217662e-19;
        out.FACstd(tt+istep) = sqrt(out.dE(Eind_fac) * out.IeCov(Eind_fac,Eind_fac,tt+istep) ...
            * out.dE(Eind_fac)')*1.60217662e-19;

        % Power carried by the precipitating electrons
        EdE = out.dE(Eind_fac)'.*out.Ec(Eind_fac)';
        out.Pe(tt+istep) = sum(out.Ie(Eind_fac,tt+istep).*EdE)*1.60217662e-19;
        out.PeStd(tt+istep) = sqrt( EdE' * out.IeCov(Eind_fac, ...
                                                    Eind_fac,tt+istep) ...
                                    * EdE ) * 1.60217662e-19;

        % Nflux, Eflux and E0
        out.Nflux(:,tt+istep) = out.Ie(:,tt+istep);
        out.Eflux(:,tt+istep) = out.Nflux(:,tt+istep) .* out.Ec(:);
        [EfMax,iMax] = max(out.Eflux(:,tt+istep));
        out.E0(tt+istep) = out.Ec(iMax);

    end

    tstamp = datestr(datenum(datetime((out.te(tt) + out.ts(tt))/2,'ConvertFrom','posixtime')));
    % print E0 only if there is some enough precipitation to fit a reasonable spectrum
    if out.Pe(tt)*1e3 < .1
        fprintf(['%s %i order, chisqr: %3.1f, FAC: %5.2f (%5.2f) uA/m2, Power: %5.2f (%5.2f) mW/m2\n'],tstamp,locminAIC,out.chisqr(tt),out.FAC(tt)*1e6,out.FACstd(tt)*1e6,out.Pe(tt)*1e3,out.PeStd(tt)*1e3)
    else
        fprintf(['%s %i order, chisqr: %3.1f, FAC: %5.2f (%5.2f) uA/m2, Power: %5.2f (%5.2f) mW/m2, E0: %5.2f keV\n'],tstamp,locminAIC,out.chisqr(tt),out.FAC(tt)*1e6,out.FACstd(tt)*1e6,out.Pe(tt)*1e3,out.PeStd(tt)*1e3,out.E0(tt)/1000)
    end        


    % update the output file
    % do not save the large density covariance matrices
    ElSpecOut = out;
    ElSpecOut.neEndCov = [];
    % save the large flux covariance only if the user explicitly
    % asked for that
    if ~out.saveiecov
        ElSpecOut.IeCov = [];
    end

    if tt==1 | any(mod(tt:tt+out.nstep-1,ndtsave)==0 )
        try
            save(out.outputfile,'ElSpecOut','-v7.3');
        catch
            ;
        end
    end

    % plot the data
    if out.plotres
        if tt == 1
            plotfig = ElSpecPlot(ElSpecOut);
        else
            if any(mod(tt:tt+out.nstep-1,ndtsave)==0)
                ElSpecPlot(ElSpecOut,plotfig);
            end
        end
    end
end
try
    save(out.outputfile,'ElSpecOut','-v7.3');
    efig = ElSpecPlot(ElSpecOut,'visible','off');
    [part1,part2,~] = fileparts(out.outputfile);
    figfile = fullfile(part1,part2);
    print(efig,figfile,'-dpng');
    close(efig)
catch
    ;
end
if out.plotres
    ElSpecPlot(ElSpecOut,plotfig);
end

end
