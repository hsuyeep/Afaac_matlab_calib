% Script to carry out a confusion noise analysis of images generated using
% the 40-channel (120KHz) and 1sec integrated data, from all 5 subbands.
% INCOMPLETE.. TODO
% pep/19May12

% Load a "clean" timeslice from the SB00?_1415-1420.bin files
fnames = {'../rawdat/SB000_1415-1420.bin', '../rawdat/SB001_1415-1420.bin', '../rawdat/SB002_1415-1420.bin', '../rawdat/SB003_1415-1420.bin', '../rawdat/SB004_1415-1420.bin'}; 
% fnames = {'SB000_1415-1420.bin'};
nfiles = length(fnames);
disp (['Found ' num2str(nfiles) ' files']);
acm = complex (zeros (288, 288, nfiles));
tobs = zeros (1, nfiles);
freq = zeros (1, nfiles);
clean_off = 10;

calvis = complex (zeros (288, 288, nfiles));
gainsol = complex (zeros (288, 1, nfiles));
uvflag = zeros (288, 288); % NOTE: Currently, the same flags go to all channels

duv = 2;
Nuv = 1000;                    % size of gridded visibility matrix
uvpad = 1024;                  % specifies if any padding needs to be added
calmap = zeros (uvpad, uvpad, nfiles);

% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal'); 

% Generate uv coordinates in local horizon coordinate system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
% normal to CS002 field (ITRF)
normal = [0.598753, 0.072099, 0.797682].'; 

for ind = 1:nfiles
	[acm(:,:,ind), tobs(ind), freq(ind)] = readms2float (fnames{ind}, clean_off, -1);
end


% Carry out a calibration of each timeslice
for ind = 1:nfiles
	disp(['Calibrating time: ' num2str(tobs(ind)) ' Freq: ' num2str(freq(ind))]);
	[calvis(:,:,ind), gainsol(:,ind), sigmas, sigman, good] = pelican_sunAteamsub (acm (:,:,ind), tobs(ind), freq(ind), uvflag, 2);

    % Generate per-subband images from each subband
    [calmap(:,:,ind), ~] = fft_imager (calvis(:,:,ind), uloc, vloc, duv, Nuv, uvpad);
end

% Integrate visibilities across subbands, print estimate of bandwidth 
% decorrelation expected due to integration

% Generate integrated image over 5 subbands

% Generate image by integrating over the image domain

% Generate noise image, SNR image
