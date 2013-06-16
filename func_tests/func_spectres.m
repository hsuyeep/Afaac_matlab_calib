% Script to determine the ideal spectral resolution from the calibration
% perspective. Done by comparison of the SNR of calibration solutions 
% generated on data integrated to different spectral resolutions.
% NOTE: Data is assumed generated using the MSToMatlab tool.
% NOTE: Data is calibrated using convergent calibration on a single 
%       timeslice.
% pep/11Jun13

% Subband 0
load '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB000_LBA_OUTER_allch_6sec_Rxxb0t0000.mat';
load '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB000_LBA_OUTER_allch_6sec_info.mat';
SB000_Rxxb0t0000 = Rxxb0t0000;
SB000_Frequencies = Frequencies;
SB000_TimeSlots = TimeSlots;

% For LBA_OUTER_BAND60 data
flagant = [1:12, 51, 193, 239,273 ]; 

addpath ../
clear pelican_sunAteamsub;
uvflag = eye(288);
debuglev = 0;
ptSun = 1;
posfilename = 'poslocal_outer.mat';
nfreq = length (SB000_Frequencies);

%%%%%%%%%%%%%%%%%%% Pre calibration plots %%%%%%%%%%%%%%%%%%%%
% Raw Bandpass, to check for systematics or RFI in input data.
figure;
for stat = 1:6
	subplot (2,3,stat);
	for ind = 1:48
		dip = (stat-1)*48 + ind;
		plot (squeeze (SB000_Rxxb0t0000(dip, dip, :)), '.-');
		hold on;
	end;
	title (sprintf ('CS00%d Chan0: %.2f\n', stat+1, SB000_Frequencies(1)));
	xlabel ('Channel'); ylabel ('Autocorr');
end;
saveas (gcf, 'rawautocorr.png', 'png');

%%%%%%%%%%%%%%%%%%% Calibration %%%%%%%%%%%%%%%%%%%%
% Calibrate every channel at 3KHz resolution
for ind = 2:nfreq
	sol3KHz(ind) = pelican_sunAteamsub (SB000_Rxxb0t0000(:,:,ind), SB000_TimeSlots, SB000_Frequencies(ind), uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
end;

% Calibrate with 6KHz resolution
k=1;
for ind = 1:2:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:1 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, 'ind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol6KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

% Calibrate with 12KHz resolution
k=1;
for ind = 1:4:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:3 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, '\nind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol12KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

% Calibrate with 24KHz resolution
k=1;
for ind = 1:8:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:7 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, '\nind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol24KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

% Calibrate with 48KHz resolution
k=1;
for ind = 1:16:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:15 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, '\nind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol48KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

% Calibrate with 96KHz resolution
k=1;
for ind = 1:32:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:31 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, '\nind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol96KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

% Calibrate with 192 KHz resolution
k=1;
for ind = 1:64:nfreq
	acc = SB000_Rxxb0t0000(:,:,ind); 
	freq = SB000_Frequencies (ind);
	for facc = 1:63 % Average 2 channels
		acc = (acc + SB000_Rxxb0t0000(:,:,ind+facc))/2;
		freq = (freq + SB000_Frequencies (ind+facc))/2;
	end;
	fprintf (1, '\nind = %d, freq = %.2f, k = %d\n', ind, freq, k);
	sol192KHz(k) = pelican_sunAteamsub (acc, SB000_TimeSlots, freq, uvflag, ... 
					flagant, debuglev, ptSun, [], [], posfilename);
	k = k+1;
end;

%%%%%%%%%%%%%%%%%%% Post calibration plots %%%%%%%%%%%%%%%%%%%%
% Generate a mean solution over full subband per antenna 
gains = sol3KHz(1).gainsol;
for ind = 2:length (sol3KHz)
	gains = gains + sol3KHz(ind).gainsol;
end;
gains = gains / length (sol3KHz);
sbgain(:,1) = gains;

gains = sol6KHz(1).gainsol;
for ind = 2:length (sol6KHz)
	gains = gains + sol6KHz(ind).gainsol;
end;
gains = gains / length (sol6KHz);
sbgain(:,2) = gains;

gains = sol12KHz(1).gainsol;
for ind = 2:length (sol12KHz)
	gains = gains + sol12KHz(ind).gainsol;
end;
gains = gains / length (sol12KHz);
sbgain(:,3) = gains;

gains = sol24KHz(1).gainsol;
for ind = 2:length (sol24KHz)
	gains = gains + sol24KHz(ind).gainsol;
end;
gains = gains / length (sol24KHz);
sbgain(:,4) = gains;

gains = sol48KHz(1).gainsol;
for ind = 2:length (sol48KHz)
	gains = gains + sol48KHz(ind).gainsol;
end;
gains = gains / length (sol48KHz);
sbgain(:,5) = gains;

gains = sol96KHz(1).gainsol;
for ind = 2:length (sol96KHz)
	gains = gains + sol96KHz(ind).gainsol;
end;
gains = gains / length (sol96KHz);
sbgain(:,6) = gains;

gains = sol192KHz(1).gainsol;
for ind = 2:length (sol192KHz)
	gains = gains + sol192KHz(ind).gainsol;
end;
gains = gains / length (sol192KHz);
sbgain(:,7) = gains;

col = {'b-', 'm-', 'r-', 'k-', 'g-', 'y-', 'c-', 'w-'};
% Post correlation and calibration averaging of solutions.
% Plot mean calibration solutions over full integrated subband, offset them for better clarity
for ind = 1:7
	plot (abs(sbgain (:,ind))+ind, char(col(ind)));
	hold on;
end;
title (sprintf ('Post calibrated average solutions over 1 subband'));
xlabel ('Ant'); ylabel ('Amp');
legend ('3KHz', '6KHz', '12KHz', '24KHz', '48KHz', '96KHz', '192KHz');
saveas (gcf, 'postcal_avg_amp.png', 'png');

% Plot mean solution phases.
figure;
for ind = 1:7
	plot (angle(sbgain (:,ind))+0.2*(ind-1), char(col(ind)));
	hold on;
end;
title (sprintf ('Post calibrated average solutions over 1 subband'));
xlabel ('Ant'); ylabel ('Ph');
legend ('3KHz', '6KHz', '12KHz', '24KHz', '48KHz', '96KHz', '192KHz');
saveas (gcf, 'postcal_avg_ph.png', 'png');

%%%%%%%%%%%%%%%%%%% Calibration Residuals & iterations %%%%%%%%%%%%%%%%%%%%
% Plot calibration residuals as a function of spectral bandwidth 
% Major cycle residuals
figure;
for ind = 1:length (sol3KHz)
	plot (ind, sol3KHz(ind).pinv_sol, 'co');
	hold on;
end;
for ind = 1:length (sol6KHz)
	plot (2*ind, sol6KHz(ind).pinv_sol, 'r.');
	hold on;
end;
for ind = 1:length (sol12KHz)
	plot (4*ind, sol12KHz(ind).pinv_sol, 'c+');
	hold on;
end;
for ind = 1:length (sol24KHz)
	h = plot (8*ind, sol24KHz(ind).pinv_sol, 'g*');
	set (h, 'markersize', 12);
	hold on;
end;
for ind = 1:length (sol48KHz)
	h = plot (16*ind, sol48KHz(ind).pinv_sol, 'md');
	set (h, 'markersize', 13);
	hold on;
end;
for ind = 1:length (sol96KHz)
	h = plot (16*ind, sol96KHz(ind).pinv_sol, 'kp');
	set (h, 'markersize', 15);
	hold on;
end;
xlabel ('Channel'); ylabel ('Major cycle residual');
saveas (gcf, 'Majorcycle_all_resid.png', 'png');

figure;
pinv = [sol3KHz(32).pinv_sol sol6KHz(16).pinv_sol sol12KHz(8).pinv_sol sol24KHz(4).pinv_sol sol48KHz(2).pinv_sol sol96KHz(1).pinv_sol];
chanwidth = [3 6 12 24 48 96];
plot (chanwidth, pinv, 'r*-');
xlabel ('channel width (KHz)'); ylabel ('Major cycle residual');
saveas (gcf, 'Majorcycle_resid.png', 'png');

iter = [sol3KHz(32).calext_iters sol6KHz(16).calext_iters sol12KHz(8).calext_iters sol24KHz(4).calext_iters sol48KHz(2).calext_iters sol96KHz(1).calext_iters];
plot (chanwidth, iter, 'r*-');
xlabel ('channel width (KHz)'); ylabel ('Major cycle iterations');
saveas (gcf, 'Majorcycle_iter.png', 'png');

% Minor cycle stuff
figure;
for ind = 1:length (sol3KHz)
	plot (ind, sol3KHz(ind).stefsol.tol, 'co');
	hold on;
end;
for ind = 1:length (sol6KHz)
	plot (2*ind, sol6KHz(ind).stefsol.tol, 'r.');
	hold on;
end;
for ind = 1:length (sol12KHz)
	plot (4*ind, sol12KHz(ind).stefsol.tol, 'c+');
	hold on;
end;
for ind = 1:length (sol24KHz)
	h = plot (8*ind, sol24KHz(ind).stefsol.tol, 'g*');
	set (h, 'markersize', 12);
	hold on;
end;
for ind = 1:length (sol48KHz)
	h = plot (16*ind, sol48KHz(ind).stefsol.tol, 'md');
	set (h, 'markersize', 13);
	hold on;
end;
for ind = 1:length (sol96KHz)
	h = plot (16*ind, sol96KHz(ind).stefsol.tol, 'kp');
	set (h, 'markersize', 15);
	hold on;
end;
xlabel ('Channel'); ylabel ('Minor cycle residual');
saveas (gcf, 'Minorcycle_all_resid.png', 'png');

figure;
tol = [sol3KHz(32).stefcal_tol sol6KHz(16).stefcal_tol sol12KHz(8).stefcal_tol sol24KHz(4).stefcal_tol sol48KHz(2).stefcal_tol sol96KHz(1).stefcal_tol];
plot (chanwidth, tol, 'r*-');
xlabel ('channel width (KHz)'); ylabel ('Minor cycle residual');
saveas (gcf, 'Minorcycle_resid.png', 'png');

stefiter = [sol3KHz(32).stefcal_iters sol6KHz(16).stefcal_iters sol12KHz(8).stefcal_iters sol24KHz(4).stefcal_iters sol48KHz(2).stefcal_iters sol96KHz(1).stefcal_iters];
plot (chanwidth, stefiter, 'r*-');
xlabel ('channel width (KHz)'); ylabel ('Minor cycle iterations');
saveas (gcf, 'Minorcycle_iter.png', 'png');

% Plot solutions as a function of bandpass.
for ind = 1:length (sol3KHz) gain3khz (:, ind) = sol3KHz(ind).gainsol; end;
for ind = 1:length (sol6KHz) gain6khz (:, ind) = sol6KHz(ind).gainsol; end;
for ind = 1:length (sol12KHz) gain12khz (:, ind) = sol12KHz(ind).gainsol; end;
for ind = 1:length (sol24KHz) gain24khz (:, ind) = sol24KHz(ind).gainsol; end;
for ind = 1:length (sol48KHz) gain48khz (:, ind) = sol48KHz(ind).gainsol; end;
for ind = 1:length (sol96KHz) gain96khz (:, ind) = sol96KHz(ind).gainsol; end;

% Plot spectral phase response of gain solution of a randomly choose antenna.
ant = 220;
plot (3 *[1:64], angle( gain3khz(ant,:)), 'co'); hold on;
plot (6 *[1:32], angle( gain6khz(ant,:)), 'r.'); hold on;
plot (12*[1:16], angle( gain12khz(ant,:)), 'c+'); hold on;
h = plot (24*[1: 8], angle( gain24khz(ant,:)), 'g*'); hold on;
set (h, 'markersize', 12);
h = plot (48*[1: 4], angle( gain48khz(ant,:)), 'md'); hold on;
set (h, 'markersize', 13);
h = plot (96*[1: 2], angle( gain96khz(ant,:)), 'kp'); hold on;
set (h, 'markersize', 15);
xlabel ('Frequency (KHz)'); ylabel ('Gain phase (rad)');
saveas (gcf, 'Ant_220_calibrated_gain_ph.png', 'png');


% Plot examinevis plot for uncalibrated 3KHz visibility.
load ('poslocal_outer.mat');
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
wloc = meshgrid (poslocal(:,3)) - meshgrid (poslocal (:,3)).';
[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 
uvw = [uloc_flag(:), vloc_flag(:), wloc_flag(:)];
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
visamphithresh= 1.5;% Reject visibilities with median >visampthresh*median.
visamplothresh= 0.5;% Reject visibilities with median >visampthresh*median.
% Arbit. choosing channel 32
% Show uncalibrated visbilities.
chan = 32;
acc = SB000_Rxxb0t0000 (:,:, chan);
[uvflag, flagant] = flagdeadcorr (acc, SB000_TimeSlots(1), ... 
				SB000_Frequencies (chan), visamphithresh, visamplothresh);


acc = acc .* (1-uvflag) - diag(diag(acc));
subplot (1,2,1);
plot (real (acc(:)), imag(acc(:)), '.'); xlabel ('Real'); ylabel ('Imag');
subplot (1,2,2);
plot (uvdist(:), abs (acc(:)), '.');
xlabel ('UV Dist. (m)'); ylabel ('Visibility amp');
title ('Uncalibrated visibilities: Chan 32/3KHz');
% Show calibrated visibilities for the same channel.
% TODO: NOT WORKING!!! GENERATE FLAGGED UVDIST!
figure;
acccal = sol3KHz (chan).calvis;
plot (real (acccal(:)), imag(acccal(:)), '.'); xlabel ('Real'); ylabel ('Imag');
subplot (1,2,2);
plot (uvdist(:), abs (acccal(:)), '.');
xlabel ('UV Dist. (m)'); ylabel ('Visibility amp');
title ('Calibrated visibilities: Chan 32/3KHz');


% Plot spectral phase response of a randomly chosen baseline.
a0 = 25; a1 = 220;
uncal_125_220 = SB000_Rxxb0t0000(a0, a1, :); % Grab all channels
for ind = 1:length(sol3KHz)
	cal_125_220 (ind) = sol3KHz(ind).calvis(a0, a1); % Calibrated visiblity
end;
plot (squeeze (angle (uncal_125_220)), 'r.-'); hold on;
plot (angle (cal_ph), 'b*-');
xlabel ('channel'); ylabel ('Phase (rad)');
title ('Spectral phase response of baseline 125:220');
legend ('uncalib', 'calib');
saveas (gcf, 'bline125__220_spec_phase_resp.png', 'png');

% Plot sigmas over bandpass 
figure;
for ind = 1:length (sol3KHz)
	plot (ind, sol3KHz(ind).sigmas, 'co');
	hold on;
end;
for ind = 1:length (sol6KHz)
	plot (2*ind, sol6KHz(ind).sigmas, 'r.');
	hold on;
end;
for ind = 1:length (sol12KHz)
	plot (4*ind, sol12KHz(ind).sigmas, 'c+');
	hold on;
end;
for ind = 1:length (sol24KHz)
	h = plot (8*ind, sol24KHz(ind).sigmas, 'g*');
	set (h, 'markersize', 12);
	hold on;
end;
for ind = 1:length (sol48KHz)
	h = plot (16*ind, sol48KHz(ind).sigmas, 'md');
	set (h, 'markersize', 13);
	hold on;
end;
for ind = 1:length (sol96KHz)
	h = plot (16*ind, sol96KHz(ind).sigmas, 'kp');
	set (h, 'markersize', 15);
	hold on;
end;
xlabel ('Channel'); ylabel ('CygA Flux estimate');
saveas (gcf, 'Sigmas_all.png', 'png');

figure;
sigmas = [sol3KHz(32).sigmas sol6KHz(16).sigmas sol12KHz(8).sigmas sol24KHz(4).sigmas sol48KHz(2).sigmas sol96KHz(1).sigmas];
plot (chanwidth, sigmas(2,:), 'r*-');
xlabel ('channel width (KHz)'); ylabel ('CygA Flux estimate');
saveas (gcf, 'Sigmas_ch32.png', 'png');

