function calibrateData(config)

% calibrateData(config)
%
% This function performs the automatic reduction of a complete LOFAR station
% calibration observation, including RFI and bad element detection and post-
% calibration flagging of erroneous results. The results are stored in a
% struct named "calres" with the following fields:
%
%   calx0   : Nrun x Nsb x Nelem matrix with calibration factors per run, per
%             subband and per element for the array of x-dipoles before post-
%             calibration flagging
%   caly0   : the corresponding result for the array of y-dipoles
%   calx    : Nrun x Nsb x Nelem matrix with calibration factors per run, per
%             subband and per element for the array of x-dipoles after post-
%             calibration flagging
%   caly    : the corresponding result for the array of y-dipoles
%   sigmax  : Nrun x Nsb x Nsrc matrix containing the estimated calibration
%             source powers per run, per subband and per calibration source for
%             the array of x-dipoles. The first source is used as flux
%             reference and has unit flux by definition.
%   sigmay  : the corresponding result for the array of y-dipoles
%   AoverTx : Nrun x Nsb matrix with an estimate of A/T towards Cas A for the
%             array of x-dipoles.
%   AoverTy : the corresponding result for the array of y-dipoles
%   t_obs   : Nrun x Nsb matrix with reconstructed observing times for each
%             snapshot observation
%   freq    : Nrun x Nsb matrix with center frequencies for each snapshot
%             observation
%
% This struct is used by the makeCalTable script to produce the static
% station calibration tables. For a datailed overview of the station
% calibration pipeline performed by this function, please see [1].
%
% The input argument to this function is a struct "config" with the following
% fields:
%
%   dirname   : name of the directory containing the calibration data 
%   bandstart : start frequency (in Hz) of band pass filter
%   bandstop  : stop frequency (in Hz) of band pass filter
%   startfreq : start frequency (in Hz) of the Nyquist zone
%   stopfreq  : stop frequency (in Hz) of the Nyquist zone
%   srcsel    : list of sources to be included in the source model
%   posfile   : AntennaArrays.conf file to be used
%   rcumode   : RCU mode number
%   EU        : boolean value, true for European station
%   uvmask    : mask to manually flag specific correlations
%   RFItol    : tolerance used in RFI detection
%   outfile   : name of the output file to store the results
%
% SJW, December 2010
% modified on 18 May 2011 by SJW to use ITRF coordinates
%
% References
% [1] Stefan J. Wijnholds, "Overview of the Initial Production Version of the
%     Station Calibration Pipeline", Internal report LOFAR-ASTRON-MEM-263,
%     ASTRON, Dwingeloo (The Netherlands), January 2011

restriction = 4;     % baseline restriction in wavelengths
maxrestriction = 30; % maximum baseline restriction in meters
Nfreq = 512;         % number of subbands
substation = 0;      % 0: full station, 1: substation 0, 2: substation 1

% extract parameters from input struct
Nsrc = length(config.srcsel);
startfreq = config.startfreq;
stopfreq = config.stopfreq;
Ns = (stopfreq - startfreq) / 512;

[~, lsout] = dos(['ls -1 ' config.dirname '*acc*']);
datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});

[pos, ~, normal] = parseAntennaFieldITRF(config.posfile, config.rcumode, config.EU);
Nelem = size(pos, 1) - (substation > 0) * size(pos, 1) * 0.5;
elemidxx = 1:2:2*Nelem;
elemidxy = 2:2:2*Nelem;

calres.calx0 = zeros(nfiles, Nfreq, Nelem);
calres.caly0 = zeros(nfiles, Nfreq, Nelem);
calres.calx = zeros(nfiles, Nfreq, Nelem);
calres.caly = zeros(nfiles, Nfreq, Nelem);
calres.sigmax = zeros(nfiles, Nfreq, Nsrc);
calres.sigmay = zeros(nfiles, Nfreq, Nsrc);
calres.AoverTx = zeros(nfiles, Nfreq);
calres.AoverTy = zeros(nfiles, Nfreq);
calres.t_obs = zeros(nfiles, Nfreq);
calres.freq = zeros(nfiles, Nfreq);

for idx = 1:nfiles
    disp(['Processing file number ' num2str(idx) ' of ' num2str(nfiles)]);
    fid = fopen(datafiles{1}{idx}, 'r');
    data = fread(fid, 'double');
    fclose(fid);
    refidx = strfind(datafiles{1}{idx}, '_acc_');
    if (substation == 0)
        acc = reshape(data(1:2:end) + 1i * data(2:2:end), [2 * Nelem, 2 * Nelem, Nfreq]);
    else
        acc = reshape(data(1:2:end) + 1i * data(2:2:end), [4 * Nelem, 4 * Nelem, Nfreq]);
        if (substation == 1)
            acc = acc(1:2*Nelem, 1:2*Nelem, :);
        else
            acc = acc(2*Nelem+1:4*Nelem, 2*Nelem+1, 4*Nelem, :);
        end
    end
    clear data;
    freq = startfreq:(stopfreq - startfreq)/Nfreq:stopfreq;
    freq = freq(1:end-1);
    t_obs = datenum(datafiles{1}{idx}(refidx-15:refidx-1), 'yyyymmdd_HHMMss') - ((Nfreq-1)/(24 * 3600):-1/(24 * 3600):0);
    % select usable subbands
    sbsel = freq > config.bandstart & freq < config.bandstop;
    [cleansbx, ~] = RFIdetection(acc(elemidxx, elemidxx, :), Ns, config.RFItol);
    [cleansby, ~] = RFIdetection(acc(elemidxy, elemidxy, :), Ns, config.RFItol);
    sbselx = sbsel & cleansbx.';
    sbsely = sbsel & cleansby.';

    % check for bad antennas
    if (sum(sbselx) == 0 || sum(sbsely) == 0)
        continue;
    end
    [selx, ~] = detectBadElem(acc(elemidxx, elemidxx, sbselx));
    [sely, ~] = detectBadElem(acc(elemidxy, elemidxy, sbsely));

    if (sum(selx) == 0 || sum(sely) == 0)
        continue;
    end
    
    % calibration
    calx = zeros(length(freq), sum(selx));
    sigmax = zeros(length(freq), Nsrc);
    Sigmanx = zeros(length(freq), Nelem, Nelem);
    uvmask = config.uvmask(selx, selx);
    [calx(sbselx, :), sigmax(sbselx, :), Sigmanx(sbselx, selx, selx)] = statcal(acc(elemidxx(selx), elemidxx(selx), sbselx), t_obs(sbselx), freq(sbselx), pos(selx, :), config.srcsel, normal, restriction, maxrestriction, uvmask);
    calx0 = checkCal(calx, selx, length(freq), Nelem);

    caly = zeros(length(freq), sum(sely));
    sigmay = zeros(length(freq), Nsrc);
    Sigmany = zeros(length(freq), Nelem, Nelem);
    uvmask = config.uvmask(sely, sely);
    [caly(sbsely, :), sigmay(sbsely, :), Sigmany(sbsely, sely, sely)] = statcal(acc(elemidxy(sely), elemidxy(sely), sbsely), t_obs(sbsely), freq(sbsely), pos(sely, :), config.srcsel, normal, restriction, maxrestriction, uvmask);
    caly0 = checkCal(caly, sely, length(freq), Nelem);

    % determine A/T
    AoverTx = findAoverT(calx, Sigmanx, freq.', t_obs.');
    AoverTy = findAoverT(caly, Sigmany, freq.', t_obs.');

    % collect result
    calres.calx0(idx, :, selx) = calx;
    calres.caly0(idx, :, sely) = caly;
    calres.calx(idx, :, :) = calx0;
    calres.caly(idx, :, :) = caly0;
    calres.sigmax(idx, :, :) = sigmax;
    calres.sigmay(idx, :, :) = sigmay;
    calres.AoverTx(idx, :) = AoverTx;
    calres.AoverTy(idx, :) = AoverTy;
    calres.t_obs(idx, :) = t_obs;
    calres.freq(idx, :) = freq;    
    save(config.outfile, 'calres')
end
