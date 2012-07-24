% Demonstration of AARTFAAC imaging, with calibration including removal of
% A-team sources (but not the sun).
%
% SJW/pep, 21 January 2012


%%  Name of directory holding snapshot data matrices
dirname = 'sterp_tseries/';
% [~, lsout] = dos(['ls -1 ' dirname '*clean.txt']);
[~, lsout] = dos(['ls -1 ' dirname '*clean.mat']);

datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
Nelem = 288; % NOTE Hardcoded

% Needed do be done only once: load ITRF positions of antennas
fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
posITRF = zeros(Nelem, 3);
posITRF(:, 1) = posfiledata{3};
posITRF(:, 2) = posfiledata{4};
posITRF(:, 3) = posfiledata{5};

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF * rotmat;

% provide meta-data
l = -1:0.005:1;                            % (l,m)-grid for sky images
m = -1:0.005:1;
c = 2.99792e8;                             % spped of light in m/s
lon = 6.869837540;                         % longitude of CS002 in degrees
lat = 52.915122495;                        % latitude of CS002 in degrees
freq = 59672546;                           % frequency in Hz
% tobs = datenum([2011, 9, 21, f_tobs(1), f_tobs(2), f_tobs(3)]);  % time in UTC
normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)

% specify calibration sources (Cas A, Cyg A, Tau A, Vir A and Sun)
srcsel =  [324, 283, 88, 179, 0];


%% Operate on all available files
for fnum = 2:nfiles
    
    % read in the cross correlation matrix, while saving it as a .mat file
    disp(['Processing file number ' num2str(fnum) ' of ' num2str(nfiles) '(' datafiles{1}{fnum} ')']);
    datestr (now)
    %acm = fill_acm_sterp (datafiles {1}{file}, 288);


    % Extract out the timestamp of the observation and load cross correlation
    % matrix
    %ccm_fname = 'sterp_tseries/SB000_115908_sterp_clean.mat';
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);

    [tok ccm_fname] = strtok (ccm_fname, '/');
    a = size (ccm_fname);
    while a(2) ~= 0
        [tok ccm_fname] = strtok (ccm_fname, '/');
        a = size (ccm_fname);
    end
    
    f_tobs = sscanf (tok, 'SB000_%2d%2d%2d')

% load AARTFAAC covariane matrix
%load acm_sterp
%load sterp_tseries/SB000_120908_try_clean
% the array covariance matrix is stored in upper triangular form, so we
% complete it
acc = acm + acm' - diag(diag(acm));
Nelem = size(acc, 1);

%% NOTE: Moved outside the loop
% % load ITRF positions of antennas
% fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
% posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
% posITRF = zeros(Nelem, 3);
% posITRF(:, 1) = posfiledata{3};
% posITRF(:, 2) = posfiledata{4};
% posITRF(:, 3) = posfiledata{5};
% 
% % convert positions to local horizon frame
% % rotation matrix taken from AntennaField.conf file from CS002
% rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
%            0.9928230000, -0.0954190000, 0.0720990000; ...
%            0.0000330000,  0.6030780000, 0.7976820000];
% poslocal = posITRF * rotmat;
% 
% % provide meta-data
% l = -1:0.005:1;                            % (l,m)-grid for sky images
% m = -1:0.005:1;
% c = 2.99792e8;                             % spped of light in m/s
% lon = 6.869837540;                         % longitude of CS002 in degrees
% lat = 52.915122495;                        % latitude of CS002 in degrees
% freq = 59672546;                           % frequency in Hz
tobs = datenum([2011, 9, 21, f_tobs(1), f_tobs(2), f_tobs(3)]);  % time in UTC
% normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)
% 
% % specify calibration sources (Cas A, Cyg A, Tau A, Vir A and Sun)
% srcsel =  [324, 283, 88, 179, 0];

%% calibrate the data
% whitening of the array covariance matrix (required for DOA estimation)
disp (['  --> Carrying out calibration with model source positions...']);
acc = acc ./ sqrt(diag(acc) * diag(acc).');
% calibrate using source positions from 3C catalog and a baseline
% restriction of 10 wavelengths or 60 meters (whichever is smaller)
[cal1, sigmas1, Sigman1] = statcal(acc, tobs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem));
% calibrate visibilities based on first calibration iteration.
acccal1 = (cal1' * cal1) .* (acc - squeeze(Sigman1));

% all sky image after first calibration iteration.
%skymapcal1 = acm2skyimage(acccal1, poslocal(:, 1), poslocal(:, 2), freq, l, m);

%% Actual source positions may differ from the known positions due to
% ionospheric image distortions, so we use a high-resolution DOA estimation
% technique (weighted subspace fitting, WSF) to estimate source positions
% relative to the position of Cas A.
% We only estimate source positions for sources with an apparent flux
% larger than 1% of the apparent flux of Cas A.
disp (['  --> Estimating source positions from data...']);
sel = sigmas1 > 0.01;
nsrc = sum(sel);
% use source positions from 3C catalog as initial estimate
load srclist3CR
if (sum(srcsel == 0) ~= 0)
    [raSun, decSun] = SunRaDec(JulianDay(tobs));                % J2000
    rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
    decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
else
    rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha];  % B1950
    decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta]; % B1950
end
[lsrc0, msrc0] = radectolm(rasrc(sel), decsrc(sel), JulianDay(tobs), lon, lat);

% find source positions using WSF algorithm
[v, d] = eig(acc);
[d, order] = sort(diag(abs(d)), 'descend');
v = v(:, order);
Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman1))) * eye(nsrc))^2 / diag(d(1:nsrc));
Es = v(:, 1:nsrc);
EsWEs = Es * Wopt * Es';
theta = fminsearch(@(x) WSFcostfun(x, EsWEs, diag(conj(1 ./ cal1)), freq, poslocal(:, 1).', poslocal(:, 2).'), [lsrc0; msrc0]);
lsrchat = theta(1:nsrc);
msrchat = theta(nsrc+1:end);

%% use updated source positions to improve calibration in a second
% iteration
disp (['  --> Carrying out second calibration with estimated source positions...']);
A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * lsrchat.' + poslocal(:, 2) * msrchat.'));
[ghat, sigmas2, Sigman2] = cal_ext(acc, A, sigmas1(sel), squeeze(abs(Sigman1)) > 0);
cal2 = (1 ./ ghat)';
% calibrate visibilities removing extended emission and system noise, using
% solutions derived by using estimated source positions
acccal2 = (cal2' * cal2) .* (acc - squeeze(Sigman2));

% all sky image after second calibration iteration, but before source subtraction
% disp (['  --> Generating skymap after second calibration...']);
% skymapcal2 = acm2skyimage(acccal2, poslocal(:, 1), poslocal(:, 2), freq, l, m);

%% display generated image
% mask = NaN(length(l));
% mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;
% figure
% imagesc(l, m, skymapcal .* mask);
% set(gca, 'FontSize', 16);
% title('calibrated sky map, extended emission removed');
% axis equal
% axis tight
% xlabel('South \leftarrow m \rightarrow North');
% ylabel('East \leftarrow l \rightarrow West');
% set(colorbar, 'FontSize', 16);

%% A-team subtraction
% We start by removing the A-team sources. Since the Sun is a barely
% resolved source given the AARTFAAC resolution of 0.8 degrees at 60 MHz,
% it requires separate treatment.
% The difference between the actual positions as estimated by WSF and the
% positions in the catalog is sufficiently large to make subtraction based
% on catalog positions fail, we use the source positions estimated above.
% reduce srcsel to source that appear brighter than 1% if Cas A, note that
% 0 is used to denote the Sun
disp (['  --> Carrying out A-team source subtraction with estimated source positions...']);
srcsel = srcsel(sel);
% construct model visibilities for A-team sources

% NOTE: Use this if the sun is not to be subtracted!
% A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * lsrchat(srcsel ~= 0).' + poslocal(:, 2) * msrchat(srcsel ~= 0).'));
% RAteam = A * diag(sigmas2(srcsel ~=0)) * A';

% NOTE: Use this if the sun is to be subtracted
A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * lsrchat.' + poslocal(:, 2) * msrchat.'));
RAteam = A * diag(sigmas2) * A';

% subtract A-team
accsubAteam = acccal2 - RAteam;

% all-sky map after removal of the A-team sources
% disp (['  --> Generating skymap after subtracting A-team...']);
% skymapsubAteam = acm2skyimage(accsubAteam, poslocal(:, 1), poslocal(:, 2), freq, l, m);
% figure
% imagesc(l, m, skymapsubAteam .* mask);
% set(gca, 'FontSize', 16);
% title('extended emission and A-team removed');
% axis equal
% axis tight
% xlabel('South \leftarrow m \rightarrow North');
% ylabel('East \leftarrow l \rightarrow West');
% set(colorbar, 'FontSize', 16);

%% Sun subtraction
% The Sun is modeled using a 3-sparse grid of pixels around the estimated
% position of the Sun. High resolution imaging using the MVDR beamformer to
% obtain a model of the Sun in the image domain proved unsuccessful.
% Finding the 3-sparse solution using an exhaustive search it the most
% computationally demanding part of this script. Hopefully Arash can
% improve on this using convex optimization.
% l0sun = lsrchat(srcsel == 0);
% m0sun = msrchat(srcsel == 0);
% % define grid of pixels
% res = 0.0015;
% delta = 0.006;
% [lsun, msun] = meshgrid(l0sun-delta:res:l0sun+delta, m0sun-delta:res:m0sun+delta);
% % initiallize matrix containing least squares cost value
% cost = NaN(length(lsun(:)), length(lsun(:)), length(lsun(:)));
% % pre-calculation, system noise added to ensure the covariance matrix is PD
% invR = inv(accsubAteam + (cal2' * cal2) .* diag(diag(Sigman2)));
% % Perform exhaustive search over all possible combinations of three
% % components on the grid defined above. Deconvolution is done using LS
% % imaging.
% for idx1 = 1:length(lsun(:))-2
%     disp(['current value of idx1: ' num2str(idx1)]);
%     pause(0.01); % to ensure progress indication is displayed
%     for idx2 = idx1+1:length(lsun(:))-1
%         for idx3 = idx2+1:length(lsun(:))
%             A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * [lsun(idx1), lsun(idx2), lsun(idx3)] + poslocal(:, 2) * [msun(idx1), msun(idx2), msun(idx3)]));
%             % calculate deconvolution matrix
%             M = abs(A' * invR * A).^2;
%             % calculate dirty image
%             sd = khatrirao(conj(invR * A), invR * A)' * accsubAteam(:);
%             % estimate flux
%             flux = real(M \ sd);
%             % if flux values are physically meaningful, assess quality of
%             % fit using least squares cost function
%             if sum(flux < 0) == 0
%                 cost(idx1, idx2, idx3) = sum(sum(abs(accsubAteam - A * diag(flux) * A').^2));
%             end
%         end
%     end
% end
% % find the three components giving lowest least squares cost
% [~, minidx] = min(testval(:));
% [idx1, idx2, idx3] = ind2sub([length(lsun(:)), length(lsun(:)), length(lsun(:))], minidx);
% % construct source model using these compoents
% A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * [lsun(idx1), lsun(idx2), lsun(idx3)] + poslocal(:, 2) * [msun(idx1), msun(idx2), msun(idx3)]));
% M = abs(A' * invR * A).^2;
% sd = khatrirao(conj(invR * A), invR * A)' * accsubAteam(:);
% flux = M \ sd;
% RSun = A * diag(flux) * A';
% % subtract A-team and Sun
% accsubSun = acccal - RAteam - RSun;

% all-sky map after removal of the A-team sources and the Sun
% skymapsubSun = acm2skyimage(accsubSun, poslocal(:, 1), poslocal(:, 2), freq, l, m);
% figure
% imagesc(l, m, skymapsubSun .* mask);
% set(gca, 'FontSize', 16);
% title('extended emission, Sun and A-team removed');
% axis equal
% axis tight
% xlabel('South \leftarrow m \rightarrow North');
% ylabel('East \leftarrow l \rightarrow West');
% set(colorbar, 'FontSize', 16);

%% Save important variables in the various stages of processing
outfilename = strrep (tok, '.mat', '_var.mat');
disp (['  --> Saving parameters to file: ' outfilename]);
% save (outfilename, 'cal1', 'sigmas1', 'Sigman1', 'acccal1', 'lsrchat', 'msrchat', 'cal2', 'sigmas2', 'Sigman2', 'acccal2', 'accsubAteam', 'tobs', 'freq', 'l', 'm', 'skymapcal1', 'skymapcal2', 'skymapsubAteam');
% save (outfilename, 'cal1', 'sigmas1', 'Sigman1', 'acccal1', 'lsrchat', 'msrchat', 'cal2', 'sigmas2', 'Sigman2', 'acccal2', 'accsubAteam', 'tobs', 'freq', 'l', 'm', 'skymapcal2', 'skymapsubAteam');
save (outfilename, 'cal1', 'sigmas1', 'Sigman1', 'acccal1', 'lsrchat', 'msrchat', 'cal2', 'sigmas2', 'Sigman2', 'acccal2', 'accsubAteam', 'tobs', 'freq', 'l', 'm');

    %% start next iteration with a clean workspace
    clear cal1 sigmas1 Sigman1 acccal1 lsrchat msrchat cal2 sigmas2 Sigman2 acccal2 accsubAteam tobs;
    %datestr (now)
    %close all

end % end the for loop over files