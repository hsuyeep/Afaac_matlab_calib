% Script to analyse the time series generated by func_tspectres.m
% Carries out de-trending and variance analysis of effect of spectral
% averaging on calibration solutions.
% pep/28Jun13

polydeg = 3;  % Detrending fit polynomial degree.
blksize = 10; % blocks of 10 secs.

% fnames = {'sol192KHz.mat', 'sol96KHz_2ch.mat', 'sol48KHz_4ch.mat', 'sol24KHz_8ch.mat', 'sol12KHz_16ch.mat'};

fnames = {'sol12KHz_16ch.mat', 'sol96KHz_2ch.mat'};
ampmean = zeros (length (fnames), 288);
ampstd  = zeros (length (fnames), 288);

for fil = 1:length (fnames);
	% Get data
	fprintf (1, 'Now working on file : %s\n', fnames {fil});
	load (fnames {fil});
	eval ([sprintf('sol = %s; clear %s;', strtok(fnames{fil},'.'), ...
										  strtok(fnames{fil}, '.'))]);
	tslots = size (sol, 1);
	nchan  = size (sol, 2);
	fprintf (1, 'From file: %s, Found %d channels, %d timeslices\n', ... 
			 fnames {fil}, nchan, tslots);

	ampts_smooth = zeros (288, tslots); % Per antenna smoothened timeseries.
	gainsol = zeros (288, blksize);
	
	% For each channel, divide the data into time blocks, find the var of each 
	% antenna's gain solution over the de-trended timeseries.
	% NOTE: Not doing for all channels, just choose the central channel.
	% for ch = 1:nchan
	if (nchan == 1) ch = 1; 
	else ch = nchan/2; end;

%%%% Test: Check std. of real/imag component of calib. solutions over time
%%%% NOTE: re. component needs de-trending.
 	for tblk = 1:round (tslots/blksize)-1; 
		indbeg = (tblk-1)*blksize + 1;	
		x = [indbeg:indbeg+blksize-1]; % Generate timeseries for fitting

		% Extract out the gainsolutions for each timeslice.
		% NOTE: 'Scalar index required for this type of multi-level indexing.' 
		% errors if for loop is avoided.
		for ts = 1:blksize
			gainsol (:, ts) = sol ((tblk-1)*blksize + ts, ch).gainsol;
		end;

		% Fit and subtract polynomial from calibration solutions
		for ant = 1:288
			ampts = real (gainsol(ant, :)); % For ALL antennas
			coe = polyfit (x, ampts, polydeg); 	
			ampts_smooth (ant, x) = ampts - polyval (coe, x);
		end;

		restd (:, tblk) = std (ampts_smooth (:, x), 0, 2);
		imstd (:, tblk) = std (imag (gainsol(:,:)), 0, 2);
	end;
	figure;
	plot (restd (32, :), 'b.-');
	hold on;
	plot (imstd (32, :), 'r.-');

%%%% Test: Check for stability of gain solutions
%%%% Code to look at the std. in phase solutions over time, as these do not show
%%%% any trend.
	% Simplified test for calibration stability
%	for tblk = 1:tslots
%		gainsol (:, tblk) = sol(tblk, ch).gainsol;
%	end;
%	ph_std (:, fil) = std (angle (gainsol.'));
	% re_std (:, fil) = std (real(gainsol'));
	% im_std (:, fil) = std (real(gainsol'));

%%%% Test: Check for stability of gain solutions
%%%% Code to subdivide timerange into blocks, over which a 3-deg. polynomial is
%%%% fitted, before subtracting the fit for detrending
	% Currently ignoring the last blksize -1 timeslots
%	for tblk = 1:round (tslots/blksize)-1; 
%		indbeg = (tblk-1)*blksize + 1;	
%		x = [indbeg:indbeg+blksize-1]; % Generate timeseries for fitting
%
%		% Extract out the gainsolutions for each timeslice.
%		% NOTE: 'Scalar index required for this type of multi-level indexing.' 
%		% errors if for loop is avoided.
%		for ts = 1:blksize
%			gainsol (:, ts) = sol ((tblk-1)*blksize + ts, ch).gainsol;
%		end;
%
%		% Fit and subtract polynomial from calibration solutions
%		for ant = 1:288
%			ampts = abs (gainsol(ant, :)); % For ALL antennas
%			coe = polyfit (x, ampts, polydeg); 	
%			ampts_smooth (ch, ant, x) = ampts - polyval (coe, x);
%		end;
%	end;
%
%	% Find variance of calibration solutions for each antenna, for this channel.
%	for ant = 1:288
%		ampmean(fil, ant)= mean(ampts_smooth (ch, ant, :)); % Mean calsolu. amp.
%		ampstd(fil, ant)= std(ampts_smooth (ch, ant, :)); % std. of calsol amp.
%		% Phase??
%	end;
end;

% Plot all errors
%col = {'b-', 'm-', 'r-', 'k-', 'g-', 'y-', 'c-', 'w-'};
%for fil = 1:length (fnames)
%	figure;
%	plot (ampstd (fil, :), char(col(fil)));
%	title (sprintf ('Gain solution amplitude error for %s', fnames{fil}));
%	xlabel ('Antenna number'); 
%	ylabel (sprintf ('std. over %d timeslices', tslots));
%	
%end;
%
