% Script for analysis of timeseries of snapshot images from AARTFAAC
% with tracking calibration.
% pep/08Feb12


% local settings, these don't change over time
freq = 59951782.2;            % TODO: SHOULD COME FROM THE BINARY DATA!
C = 299792458; % speed of light
l = [-1:0.005:1]; m = l;
restriction = 20;             % Avoid uv values with <10 wavelengths length
maxrestriction = 60;          % Avoid uv values with >60 wavelengths length
maxiter = 1;
window = 5;                  % Number of timeslices we will hold in memory

% FFT Imaging related
duv = 600/511;                % In meters, grid spacing in fft imager
Nuv = 512;                    % size of gridded visibility matrix
uvpad = 512;                  % specifies if any padding needs to be added
calmap = zeros (Nuv, Nuv, window);
calvis = calmap;

% Get a list of files containing the data cubes
% The files are usually generated from a binary dump of relevant dataslices
% from an MS. They are generated using mstofloat.py, and containmultiple
% subbands and timeslices for a certain timerange.
% NOTE: These should come from some kind of parameter text file!
fnames = {'SB000_1415-1420.bin', '1222_1227.bin'};
nfiles = 2;
Nelem = 288;                 % TODO: SHOULD COME FROM BINARY DATA
acm_store = zeros (Nelem, Nelem, window);
ghat_store = zeros (Nelem, 1, window);
acm = zeros (Nelem);
nblines = Nelem * (Nelem + 1)/2; 

% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal'); 

% Generate uv coordinates in local horizon coordinate system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
% normal to CS002 field (ITRF)
normal = [0.598753, 0.072099, 0.797682].'; 
% Generate uvw coordinates in ITRF coordinates, used for calibration
u = meshgrid(posITRF(:, 1)) - meshgrid(posITRF(:, 1)).';
v = meshgrid(posITRF(:, 2)) - meshgrid(posITRF(:, 2)).';
w = meshgrid(posITRF(:, 3)) - meshgrid(posITRF(:, 3)).';
uvw = [u(:), v(:), w(:)];

% Create a visibility flag based on wavelength restriction on (u,v)
uvflag = eye (Nelem);
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
mask = reshape(uvdist, [Nelem, Nelem]) < min([restriction * (C / freq), maxrestriction]);
mask = mask | uvflag;

% specify calibration sources from the 3CR catalog (Cas A, Cyg A, Tau A,
% Vir A and Sun)
load srclist3CR
srcsel =  [324, 283, 88, 179, 0];
sigmas1_store = zeros (length (srcsel), 1, window);
%
fid = fopen (fnames{1}, 'rb');

slnum = 0;
matrixify = triu (ones (Nelem));  % To convert incoming complex vector to matrix

% % Call a full fledged iterative calibration for the first time
% % calibrate the data for the current snapshot, using no initial
% % estimates
% % Read in the timestamp for this timeslice, available as MJD in seconds
% 	% from the measurement set.
% 	t_obs1 = fread (fid, 1, 'double');
% 	disp (['Timeslice: ' num2str(t_obs1, '%f') ' (MJD secs), of ' num2str(slnum)]);
% 	t_obs = t_obs1/86400 + 2400000.5; 
% 	% Reading real and imaginary, available as a stream of floats.
% 	% even floats being real parts, odd floats being imag
% 	a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
% 	comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
% 	% create instantaneous ccm from vector
% 	acm = triu (ones (Nelem));
%     acm (acm == 1) = comp;
% 	%acm (matrixify == 1) = comp;
% 	%acm = acm + acm' - diag (diag(acm));
% 	acc = acm + acm' - diag(diag(acm));
%     
% [cal0, sigmas0, Sigman0] = statcal(acc, t_obs, freq, posITRF, srcsel, normal, 4, 20, eye(Nelem));
% disp ('Initial estimates of model src flux:'); disp (sigmas0);   
sigmas0 = [1];   
    
    %% Start the tracking part, with the initial estimates obtained
    tic
for i = [1:window] % TODO: REPLACE WITH NUMBER OF TIMESLICES
	slnum = slnum + 1;
	if (slnum >= window)
		slnum = 1;
	end

	% Read in the timestamp for this timeslice, available as MJD in seconds
	% from the measurement set.
	t_obs1 = fread (fid, 1, 'double');
	disp (['Timeslice: ' num2str(t_obs1, '%f') ' (MJD secs), of ' num2str(slnum)]);
	% Convert back to Julian day, as desired by track_statcal. NOTE:
	% JulianDay () should not be called now!
	t_obs = t_obs1/86400 + 2400000.5; 
    % t_obs_datenum = datenum([2011, 9, 21, 12, 22, 01]);

	% Reading real and imaginary, available as a stream of floats.
	% even floats being real parts, odd floats being imag
	a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
	comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
	
	%% create instantaneous ccm from vector
	acm = triu (ones (Nelem));
    acm (acm == 1) = comp;
	%acm (matrixify == 1) = comp;
	%acm = acm + acm' - diag (diag(acm));
	acc = acm + acm' - diag(diag(acm));
    
	
	%% Call the tracking calibration
	%[cal1, sigmas1, Sigman1] = tracking_cal (acc, t_obs, freq, posITRF, srcsel, srclist3CR,  normal, mask); %, cal0, sigmas0, Sigman0, uvw);
	 [ghat, sigmahat, iter, cal1, sigmas1, Sigman1, Sigmanhat] = track_statcal (acc, t_obs, freq, posITRF, srcsel, normal, 4, 20, eye(Nelem), maxiter, srclist3CR, sigmas0);
     %sigmas0 = sigmas1;
	 %disp ('Estimated model src flux: '); disp (sigmas0);
	%%
	% acm_store(:,:,slnum) = (cal1' * cal1) .* acm; % Generate calibrated visibilities 
	acccal = (cal1' * cal1) .* acm; % Generate calibrated visibilities 
	% map = acm2skyimage (acmcal, poslocal(:,1), poslocal(:,2), freq, l, m);
	% [uncalmap, uncalvis] = fft_imager_sjw (acc (:), uloc(:), vloc(:),
	% duv, Nuv, uvpad);
    ghat_store (:,:,slnum) = ghat (:, iter);
    sigmas1_store (:,slnum) = sigmas1;
	[calmap(:,:,slnum), calvis(:,:,slnum)] = fft_imager_sjw (acccal (:), uloc(:), vloc(:), duv, Nuv, uvpad);
    
end
toc