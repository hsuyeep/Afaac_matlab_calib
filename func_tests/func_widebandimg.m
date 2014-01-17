% Script to generate images by collapsing over a subband or higher 
% pep/14Jun13

%% Raw data container
acc = zeros (288, 288, 64*5); % assumed spectral resolution is 64 channels.
freq = zeros (1, 64*5);
tobs = zeros (1, 5);

% Read in all 5 subbands.
for ind = 1:5
	fprintf (2, 'Reading subband %d\n', ind-1);
	load (sprintf ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB00%d_LBA_OUTER_allch_6sec_Rxxb0t0000.mat', ind-1));
	acc(:,:,(ind-1)*64 + [1:64]) = Rxxb0t0000;
	load (sprintf ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB00%d_LBA_OUTER_allch_6sec_info.mat', ind-1));
	freq ((ind-1)*64 + [1:64]) = Frequencies;
	tobs (ind) = TimeSlots;
end;
fprintf (2, 'Timeslots read: %f \n', tobs);

% For LBA_OUTER_BAND60 data
flagant = [1:12, 51, 193, 239,273 ]; 
addpath ../
clear pelican_sunAteamsub;
uvflag = eye(288);
debuglev = 0;
ptSun = 1;
posfilename = 'poslocal_outer.mat';
% Local horizon based coordinates of array in ITRF
load (posfilename, 'posITRF', 'poslocal'); 
nfreq = length (freq)

% Imaging parameters
radec = 0;
duv = 2.5;						% Default, reassigned from freq. of obs. to
								% image just the full Fov (-1<l<1)
Nuv = 500; %1000                % size of gridded visibility matrix
uvpad = 512; %1024              % specifies if any padding needs to be added
dl = (299792458./(freq * uvpad * duv)); % dimensionless, in dir. cos. units
    
% NOTE: Total imaged Field of View is determined by the visibility 
% grid-spacing, duv.
lmax = dl * uvpad / 2;
% Local horizon based coordinates of array in ITRF
load (posfilename, 'posITRF', 'poslocal'); 

% Generate uv coordinates in local horizon coord. system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 

%% Calibrate all channels at highest resolution, skip DC of each subband
for ind = 2:nfreq
	if (mod (ind,64) == 1)
		continue;
	else
		fprintf (2, '\nCalibrating channel %d, freq %.2f\n', ind, freq(ind));
		sol3KHz(ind) = pelican_sunAteamsub (acc (:,:,ind), tobs(1), ... 
		   freq (ind), uvflag, flagant, debuglev, ptSun, [], [], posfilename);

		fprintf (2, '\nImaging channel %d, freq %.2f\n', ind, freq(ind));
		[radecskymap, lmskymap(:,:,ind), vispad, l(:,ind), m(:,ind)] = fft_imager_sjw_radec (sol3KHz(ind).calvis(:), uloc_flag(:), vloc_flag(:), duv, Nuv, uvpad, tobs(1), freq(ind), radec);
		
	end;
end;

% Generate images by filling uvplane with calibrated visibilities from all
% channels in a subband.
% Determine the uv-spacing required, based on frequency
for ind = 1:nfreq

    
end;

% Display high resolution images
for ind = 2:size (lmskymap, 3)
	imagesc (l(:,ind), m(:,ind), abs (lmskymap(:,:,ind)));
	pause (0.2);
end;

% Generate stacked image via 'simple' stacking
stackmap = zeros (512);
for ind = 2:size (lmskymap, 3)
	stackmap = stackmap + lmskymap (:,:,ind);
end;

% plot bandpass over all subbands
figure;
for stat = 1:6
	subplot (2,3,stat);
	for ind = 1:48
		dip = (stat-1)*48 + ind;
		plot (squeeze (acc (dip, dip, :)), '.-');
		hold on;
	end;
	title (sprintf ('CS00%d Chan0: %.2f\n', stat+1, freq (1)));
	xlabel ('Channel'); ylabel ('Autocorr');
end;
saveas (gcf, 'raw5sbautocorr.png', 'png');


% plot phase variation over array, for given frequency within subband
uloc_spect = zeros (length (uloc_flag), nfreq);
uvw = [uloc_flag(:), vloc_flag(:), zeros (length(uloc_flag(:)), 1)];
normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
for ind = 0:30:120
    lambda = 299792458/freq(1);
	ph = ((uvw(:,1)/lambda)*sind(ind) +  (uvw(:,2)/lambda)*sind(ind));
	% uloc_spect (:,ind) = uloc_flag (:)/freq (ind);
	plot (uvdist/lambda, ph, 'd');
    xlabel ('uvdist (lambda)'); ylabel ('Phase(rad)');
    title (sprintf ('Phase ramp across array for look angle: %d deg\n', ind));
	pause;
end;

% Plot variation of phase for a given look angle, as a function of
% frequency
for ind = 1:nfreq
    lambda = 299792458/freq(ind);
    los = 60; % deg
    sel = uvdist/lambda > 60; % Wavelength cutoff for longest blines.
	ph = ((uvw(sel,1)/lambda)*sind(los) +  (uvw(sel,2)/lambda)*sind(los));
	% uloc_spect (:,ind) = uloc_flag (:)/freq (ind);
	plot (uvdist(sel)/lambda, ph, 'd');
    xlabel ('uvdist (lambda)'); ylabel ('Phase(rad)');
    title (sprintf ('Phase ramp across array for freq: %f, look angle: %d deg\n', freq(ind), los));
	pause (0.1);
end;
