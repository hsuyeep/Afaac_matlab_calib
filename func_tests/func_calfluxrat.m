% Script to generate the time-series of the peak and integrated flux of 
% specified calibrator sources using extractcal.m. 
% Plots the peak and integ. flux timeseries, the ratios with the first calsrc,
% Also carries out outlier removal via a median filter over a specifiable 
% time window.
% pep/12Jul12

% Arguments:
%	fname : File containing calibrated images (use rdimg2bin to read in).
%  offset : Timeslice offset from which to begin.
% ntslices: Number of timeslices to process.


function func_calfluxrat (fname, offset, ntslices, posfilename, debug)

	addpath ../srcfind 				   	% Add fit2dgauss () to path.
	
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('func_calfluxrat: fid < 0! Quitting.');
		return;
	end;
		
	% Load various meta data
    rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
    rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
    rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	load (posfilename, 'posITRF', 'poslocal');
    rodata.posITRF = posITRF;
    rodata.poslocal = poslocal;
	load srclist3CR.mat

	if (ntslices < 0)
		% Figure out number of records via kludgy filesize method.
		img = readimg2bin (fid);
		fdir = dir (fname);
		filesize = fdir.bytes;
		imgwhos = whos ('img');
		imgsize = imgwhos.bytes;
		ntslices = round (filesize/imgsize);
		fprintf (1, '--> Found %d timeslices in file.\n', ntslices);
		fprintf (1, '--> NOTE: Does not work with multiple channels yet!\n');
	end;

	% NOTE: Should be an all-sky callist which is pruned to the current FoV.
	% NOTE: Current list assumes 3CR catalog!
	% calsrcs = 3C[
	callist = [200, 133, 237];
	ncal = length (callist);

	% Create datastructures
	tpk = zeros (ncal, 1); tint = zeros (ncal, 1); % Instantaneous peak/int. fl.

	peakfl = zeros (ncal, ntslices);    % Timeseries datastructures
	intfl = peakfl;
	resid = peakfl;

	winsize = 10; % Timeslices          % Window datastructures
	pkwin = zeros (ncal, winsize);
	intwin = pkwin;
	thresh = 3;                         % Rejection threshold = 10*median value
	badfit = zeros (ncal, ntslices);	% TO hold the number of badly fitted src
	badtimes = 0;						% Variable storing number of tslices 
										% discarded due to median filter.
	leg_str = cell (1, ncal);			% Cell array for dynamic legends


	% Fill in timeslice window.
	fprintf (2, '--> Filling filter window of %d timeslices...', winsize);
	for ind = 1:winsize
		img = readimg2bin (fid);
		fprintf (1, '\nRec %4d ', ind);
		[pkwin(:, ind), intwin(:, ind), upcals, resid(:, ind), exitfl] = ...
			extractcal (img.map, img.l, img.m, img.tobs, ...
					img.freq, rodata, srclist3CR, callist, true, 0);	

		badfit (:, ind) = exitfl;   % Record locations of bad fits.
		% Check for convergence of fminsearch
		% A 0 implies Max. num. of func. evals, reached. 1s is good. -1 is bad
		if (sum(exitfl) ~= ncal) 
			fprintf (1, 'Exit fl: %s\n', num2str(exitfl));
			% badfitcals = find (exitfl == 0); % Indices of inconvergent cals.
			% badfit (badfitcals) = badfit (badfitcals) + 1;
		end;
		
		% Also filling appropriate location in full timeseries
		peakfl (:, ind) = pkwin (:, ind);   intfl (:, ind) = intwin (:, ind);
		
	end;
	% Generate initial statistics for each calsource, over winsize timeslices.
	% Transpose due to median being a row vector computed over cols.
	pkmed = median (pkwin'); intmed = median (intwin'); 
	fprintf (2, '--> Median peak   flux: %s\n', num2str(pkmed, '%.2f '));
	fprintf (2, '--> Median integ. flux: %s\n', num2str(intmed, '%.2f '));
	fprintf (2, '--> Threshold: %d, %s Jy\n', thresh);

	% Move back to beginning of file.
	% TODO: But be careful about variables like badfit etc.

	for ind = winsize:ntslices
		img = readimg2bin (fid);

		fprintf (1, '\nRec %4d ', ind);
		[tpk, tint, upcals, resid(:, ind), exitfl] = ...
			extractcal (img.map, img.l, img.m, img.tobs, img.freq, rodata, ...
						 srclist3CR, callist, true, 0);	
		badfit (:, ind) = exitfl;   % Record locations of bad fits.

		% Check for convergence of fminsearch
		% A 0 implies Max. num. of func. evals, reached. 1s is good. -1 is bad
		if (sum(exitfl) ~= ncal) 
			fprintf (1, 'Exit fl: %s\n', num2str(exitfl));
			% badfitcals = find (exitfl == 0); % Indices of inconvergent cals.
			% badfit (badfitcals) = badfit (badfitcals) + 1;
		end;

		% Decide if current timeslice should enter timeseries.
		% NOTE: Not using the integrated flux!! 
		if (sum (tpk > thresh*pkmed) > 0) % Vectorial comparison.
			fprintf (2, '<-- Discarding time %.2f\n', img.tobs);	
			badtimes = badtimes + 1;
			continue;
		end;

		% Fill window circularly
		winind = mod (ind, winsize) + 1;
		pkwin (:, winind) = tpk; intwin (:, winind) = tint;
		peakfl (:, ind) = tpk;   intfl (:, ind) = tint;
	end;

	% Print statistics
	fprintf (1, '\nBad fit count: %s\n', num2str(sum (badfit' == 0), ' %2d '));
	fprintf (1, 'Discarded timeslices: %d\n', badtimes);

	% Plot fitted fluxes
	col = {'b*-', 'm*-', 'r*-', 'k*-', 'g*-', 'y*-', 'w*-', 'c*-'};
	for ind = 1:ncal
		leg_str {ind} = srclist3CR(callist(ind)).name;
	end;

	figure;
	% Peak fluxes and ratios
	subplot (221);
	for ind = 1:ncal
		% pfl = peakfl (ind, badfit(ind, :) == 0);
		% plot (pfl, char(col(ind)));
		plot (peakfl (ind, :), char(col(ind))); % Unfiltered data.
		hold on;
	end;
	title ('Extract peak flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);
	subplot (222);
	for ind = 1:ncal
		plot (peakfl(ind, :)./peakfl(1,:), char(col(ind)));
		% plot (resid(ind,:), char(col(ind)));
		hold on;
	end;
	title ('Extract peak flux ratio of calibrators'); 
	% title ('resid flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);

	% Integrated fluxes and ratios
	subplot (223);
	for ind = 1:ncal
		% ifl = intfl (ind, badfit(ind, :) == 0);
		% plot (ifl, char (col(ind)));
		plot (intfl(ind, :), char (col(ind)));
		hold on;
	end;
	title ('Extract int. flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);
	subplot (224);
	for ind = 1:ncal
		plot (intfl(ind, :)./intfl(1,:), char(col(ind)));
		hold on;
	end;
	title ('Extract int. flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);

	% Plot histograms. We discard outliers for individual timeseries histograms
	% For ratios, we discard all bad fits in the union of all cal sources. 
	figure;
	badfitunion = sum (badfit);
	pkflref = peakfl (1, badfitunion == ncal); 
	for ind = 1:ncal
		subplot (2,3,ind);
		pfl = peakfl (ind, badfit(ind, :) == 1);
		hist (pfl, 50);
		title (sprintf ('Peak flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2,3,ind+ncal);
		pfl = peakfl (ind, badfitunion == ncal); 
		hist (pfl./pkflref, 50);
		title (sprintf ('Peak flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;

	figure;
	intflref = intfl (1, badfitunion == ncal); 
	for ind = 1:ncal
		subplot (2,3,ind);
		ifl = intfl (ind, badfit(ind, :) == 1);
		hist (ifl, 50);
		title (sprintf ('Integ. flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2,3,ind+ncal);
		pfl = intfl (ind, badfitunion == ncal); 
		hist (pfl./intflref, 50);
		title (sprintf ('Peak flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;

	% Save the workspace for longer files!
	fsave = strcat (fname, '_save.mat');
	save (fsave);
