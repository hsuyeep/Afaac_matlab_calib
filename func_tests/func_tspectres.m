% Script to generate calibration solution statistics at different spectral
% resolutions. Now features error on estimates via multiple samples in time.
% NOTE that we still operate on multiple channels, a single timeslice at a time,
% but each timeslice's vital parameters are stored for generating statistics.
% pep/19Jun13

function func_tspectres ()
%% 
ntimes = 56; % Default.
nbands = 1; % Default.
debug = 0;
first = 1;

addpath ../
clear pelican_sunAteamsub;
uvflag = eye(288);
debuglev = 0;
ptSun = 1;
posfilename = 'poslocal_outer.mat';

% For LBA_OUTER_BAND60 data
flagant = [1:12, 51, 193, 239,273 ]; 


%% Read in metadata
for sb = 1:5 % NOTE: all 5 subbands being read!
	% filename for info structure of this timeslice and subband
	fname = sprintf ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB00%d/info.mat', sb-1);
	load (fname);
	if first == 1
		freq = zeros (nbands, length(Frequencies));
		tobs = zeros (nbands, ntimes);
		first = 0;
	end;
	eval ([sprintf('freq (%d, :) = Frequencies; clear Frequencies;',sb)]);
	eval ([sprintf('tim_sb%d = TimeSlots; ',sb)]);
end;
nfreq = size (freq,1) * size (freq, 2);
ntimes = length (tim_sb1); 

first = 0;

%% Raw data container
acc = zeros (288, 288, 64*nbands); % assumed spectral resolution is 64 channels,
							  % This holds data over 1 subband, 1 sec

%% Read in all timeslices
for ts = 1:ntimes
	% Read in desired subbands
	for sb = 1:nbands
		% filename for this timeslice and subband
		fname = sprintf ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/spectres/SB00%d/Rxxb0t00%02d.mat', sb-1, ts-1);
		fprintf (2, '\nReading from %s\n', fname);
		load (fname);
		eval ([sprintf('acc_sb = Rxxb0t00%02d; clear Rxxb0t00%02d;',ts-1,ts-1)]);
		acc (:,:, (sb-1)*64 + (1:64)) = acc_sb;
	end;

	% DEBUG: Show bandpass from all 5 subbands for this timeslice.
	if (debug == 1)
		for stat = 1:6
			subplot (2,3,stat);
			for ind = 1:48
				dip = (stat-1)*48 + ind;
				% plot (freq(:), squeeze(acc (dip, dip, :)), '.');
				plot (squeeze(acc (dip, dip, :)), '.');
				hold on;
			end;
			title (sprintf ('Time: %d, CS00%d\n', tim_sb1(ts), stat+1));
			xlabel ('Channel'); ylabel ('Autocorr');
		end;
	end;	

	
	% Calibrate every channel at 3KHz resolution, currently only over 1 sb.
	for ind = 1:64*nbands
		if (ind == 1) continue; end;
		sol3KHz(ts, ind) = pelican_sunAteamsub (acc (:,:,ind), TimeSlots(ts), ...
			freq(1, ind), uvflag, flagant, debuglev, ptSun, [], [], posfilename);
	end;

	fprintf (2, 'Now calibrating 6 KHz chunks...\n');
	k=1;
	for ind = 1:2:64*nbands
		acccum = acc (:,:,ind); 
		freqcum = freq (1, ind);
		for facc = 1:1 % Average 2 channels
			acccum = (acccum + acc (:,:,ind+facc))/2;
			freqcum = (freqcum + freq (1, ind+facc))/2;
		end;
 		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
		sol6KHz(ts, k) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
		k = k + 1;
	end;

%	fprintf (2, 'Now calibrating 12 KHz chunks...\n');
%	for ind = 1:4:64*nbands
%		acccum = acc (:,:,ind); 
%		freqcum = freq (1, ind);
%		for facc = 1:3 % Average 4 channels
%			acccum = (acccum + acc (:,:,ind+facc))/2;
%			freqcum = (freqcum + freq (1, ind+facc))/2;
%		end;
% 		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
%		sol12KHz(ts, ind) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
% 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
%	end;
%
%	fprintf (2, 'Now calibrating 24 KHz chunks...\n');
%	for ind = 1:8:64*nbands
%		acccum = acc (:,:,ind); 
%		freqcum = freq (1, ind);
%		for facc = 1:7 % Average 4 channels
%			acccum = (acccum + acc (:,:,ind+facc))/2;
%			freqcum = (freqcum + freq (1, ind+facc))/2;
%		end;
% 		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
%		sol24KHz(ts, ind) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
% 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
%	end;
%
%	fprintf (2, 'Now calibrating 48 KHz chunks...\n');
%	for ind = 1:16:64*nbands
%		acccum = acc (:,:,ind); 
%		freqcum = freq (1, ind);
%		for facc = 1:15 % Average 4 channels
%			acccum = (acccum + acc (:,:,ind+facc))/2;
%			freqcum = (freqcum + freq (1, ind+facc))/2;
%		end;
% 		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
%		sol48KHz(ts, ind) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
% 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
%	end;
%
%	fprintf (2, 'Now calibrating 96 KHz chunks...\n');
%	for ind = 1:32:64*nbands
%		acccum = acc (:,:,ind); 
%		freqcum = freq (1, ind);
%		for facc = 1:31 % Average 4 channels
%			acccum = (acccum + acc (:,:,ind+facc))/2;
%			freqcum = (freqcum + freq (1, ind+facc))/2;
%		end;
% 		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
%		sol96KHz(ts, ind) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
% 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
%	end;

%	for ind = 1:64:64*nbands
%		acccum = acc (:,:,ind); 
%		freqcum = freq (1, ind);
%		for facc = 1:63 % Average 4 channels
%			acccum = (acccum + acc (:,:,ind+facc))/2;
%			freqcum = (freqcum + freq (1, ind+facc))/2;
%		end;
%		fprintf (1, 'ind = %d, freq = %.2f\n', ind, freqcum);
%		sol192KHz(ts, ind) = pelican_sunAteamsub (acccum, TimeSlots(ts),...
% 		   freqcum, uvflag, flagant, debuglev, ptSun, [], [], posfilename);
%	end;

	% Generate statistics across the frequency axis on calibration solution.
%	if first == 1
%		allgain = zeros (ntimes, 64*nbands, length (sol3KHz(2).gainsol));
%		allsigmas = zeros (ntimes, 64*nbands, length (sol3KHz(2).sigmas));
%		first = 0;
%	end;
%
%	for ind = 2:length (sol3KHz)
%		allgain (ts, ind, :) = sol3KHz(ind).gainsol;
%		allsigmas (ts, ind, :) = sol3KHz(ind).sigmas;
%	end;
%
%	gains = sol3KHz(2).gainsol;
%	sigmas = sol3KHz(2).sigmas;
%	for ind = 2:64 % DOING FOR ONE SUBBAND ONLY
%		gains = gains + sol3KHz(ind).gainsol;
%		sigmas = sigmas + sol3KHz(ind).sigmas;
%	end;
%	gains = gains / length (sol3KHz);
%	sigmas = sigmas/length (sol3KHz);
%
%	% Store averaged values.
%	sbgain(:,ts) = gains;
%	sbsigmas(:,ts) = sigmas;
	
end; % End of ts loop

% save ('tspectres_1.mat', 'sol96KHz', 'sol48KHz', 'sol24KHz', 'sol12KHz', 'sol6KHz', 'sol3KHz');
save ('tspectres_2.mat', 'sol6KHz', 'sol3KHz');
