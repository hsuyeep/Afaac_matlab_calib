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
% posfilename: Array configuration of the images.
% callist : List of sources within FoV to be monitored.
% debug   : Turn on status messages.
% Returns:
%	fit   : Structure containing all relevant fit parameters.


function fit = func_calfluxrat (fname, offset, ntslices, posfilename, callist, debug)

	% addpath ../srcfind 				   	% Add fit2dgauss () to path.
	
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

	% Figure out number of records via kludgy filesize method.
	img = readimg2bin (fid);
	fdir = dir (fname);
	filesize = fdir.bytes;
	imgwhos = whos ('img');
	imgsize = imgwhos.bytes;
		
	if (ntslices < 0)
		ntslices = round (filesize/imgsize)-1; % Last timeslice is always tricky
		fprintf (1, '--> Found %d timeslices in file.\n', ntslices);
		fprintf (1, '--> NOTE: Does not work with multiple channels yet!\n');
		fseek (fid, 0, 'bof');
	end;

	% Move to desired offset
	if (offset > 1)
		% fseek (fid, offset*imgwhos.bytes, 'bof');
		fprintf (1, 'Moving to offset %d...\n', offset);
		for ind = 1:offset
			img = readimg2bin (fid);
		end;
	end;

	% NOTE: Should be an all-sky callist which is pruned to the current FoV.
	% NOTE: Current list assumes 3CR catalog!
	% calsrcs = 3C[295, 219, 338, 410.1, 10, 84, 338]
	% callist = [117, 133, 265, 256, 200, 237];  % For day time observations
	% callist = [  4, 314, 265, 256, 200, 237];  % For nite time observations
	ncal = length (callist);

	% Create datastructures
	tpk = zeros (ncal, 1); tint = zeros (ncal, 1); % Instantaneous peak/int. fl.

	tobs = zeros (ntslices,1);
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

	% Debug windows
	if (debug > 3)
		islemap = figure;
		fit2dmap = figure;
	else
		islemap = -1; fit2dmap = -1;
	end;

	% Figure out which calibraters are up via a single first call to extractcal.
	img = readimg2bin (fid);
	tobs(1) = img.tobs;
	psf = [];
	[upcals, fit(:, 1)] = ...
			extractcal (img.map, img.l, img.m, img.tobs, ...
					img.freq, rodata, srclist3CR, callist, true, psf, debug, ...
					islemap, fit2dmap);	

	% NOTE that extractcal returns an array of structures, one per fitted src.
	fitparam = zeros (ncal, length (fit(1,1).fitparams));
	% fitparam = zeros (ncal, length (firstfit(1).fitparams));
	for ind = 1:sum(upcals)
		fitparam (ind,:) = fit(ind, 1).fitparams;
		% fitparam (ind,:) = firstfit(ind).fitparams;
	end;

	fprintf (1, 'Found %d sources to monitor.', sum(upcals));

	% Create an array of structures for full timeseries.
	% fit (ncal, ntslices) = {firstfit};
	% fit (:, 1) = {firstfit};
	% TODO: Currently dynamically expanding the array of structures.


	% Fill in timeslice window.
	fprintf (2, '--> Filling filter window of %d timeslices...', winsize);
	for ind = 2:winsize
		img = readimg2bin (fid);
		tobs (ind) = img.tobs;
		fprintf (1, 'Rec %4d ', ind);
		[upcals, fit(:,ind)] = ...
			extractcal (img.map, img.l, img.m, img.tobs, ...
					img.freq, rodata, srclist3CR, callist, true, psf, debug, islemap, fit2dmap);	

		badfit (:, ind) = [fit(:,ind).exitfl];   % Record locations of bad fits.
		% Check for convergence of fminsearch
		% A 0 implies Max. num. of func. evals, reached. 1s is good. -1 is bad
		if (sum(badfit(:,ind)) ~= sum(upcals)) 
			fprintf (1, 'Exit fl: %s\n', num2str([fit(:,ind).exitfl]));
			% badfitcals = find (exitfl == 0); % Indices of inconvergent cals.
			% badfit (badfitcals) = badfit (badfitcals) + 1;
		end;
		
		% Need to use this kludge due to array of structures, each with a vector
		% member
		for cl = 1:ncal
			fitparam (cl,:) = fit(cl, ind).fitparams;
		end;
		tpk = fitparam (:,3);
		tint = 2*pi*fitparam (:,3).*fitparam (:,4).*fitparam (:,5); 

		% Also filling appropriate location in full timeseries
		pkwin  (:,ind) = tpk; 
		peakfl (:,ind) = tpk; 
		intfl  (:,ind) = tint;
		intwin (:,ind) = tint;
		tobs (ind) = img.tobs;
	end;
	% Generate initial statistics for each calsource, over winsize timeslices.
	% Transpose due to median being a row vector computed over cols.
	pkmed = median (pkwin'); intmed = median (intwin'); 
	fprintf (2, '--> Median peak   flux: %s\n', num2str(pkmed, '%.2f '));
	fprintf (2, '--> Median integ. flux: %s\n', num2str(intmed, '%.2f '));
	fprintf (2, '--> Threshold: %d, %s Jy\n', thresh);

	% Move back to beginning of file.
	% TODO: But be careful about variables like badfit etc.

	for ind = winsize+1:ntslices
		img = readimg2bin (fid);

		fprintf (1, 'Rec %4d ', ind);
		[upcals, fit(:,ind)] = ...
			extractcal (img.map, img.l, img.m, img.tobs, img.freq, rodata, ...
						 srclist3CR, callist, true, psf, debug, islemap, fit2dmap);	
		badfit (:, ind) = [fit(:, ind).exitfl];  % Record locations of bad fits.
		% Check for convergence of fminsearch
		% A 0 implies Max. num. of func. evals, reached. 1s is good. -1 is bad
		if (sum(badfit(:,ind)) ~= ncal) 
			fprintf (1, 'Exit fl: %s\n', num2str([fit(:,ind).exitfl]));
			% badfitcals = find (exitfl == 0); % Indices of inconvergent cals.
			% badfit (badfitcals) = badfit (badfitcals) + 1;
		end;

		for cl = 1:ncal
			fitparam (cl,:) = fit(cl, ind).fitparams;
		end;
		tpk = fitparam (:,3);
		tint = 2*pi*fitparam (:,3).*fitparam (:,4).*fitparam (:,5); 
		tobs (ind) = img.tobs;

		% Decide if current timeslice should enter timeseries.
		% NOTE: Not using the integrated flux!! 
		% Check for outliers and bad fits.
		badsrc = (tpk > thresh*pkmed') | (1-badfit (:,ind)); 
		% if (sum (tpk' > thresh*pkmed) > 0) % Vectorial comparison.
		if (sum (badsrc) > 0) % Vectorial comparison.
			fprintf (2, '<-- Discarding time %.2f. Badsrc: %s.\n', ...
					img.tobs, num2str(badsrc));	
			badtimes = badtimes + 1;
			tpk (badsrc == 1) = NaN; % Replace outliers with NaN.
			tint (badsrc == 1) = NaN;
			% continue;
		else

			% Fill window circularly
			winind = mod (ind, winsize) + 1;
			pkwin (:, winind) = tpk; intwin (:, winind) = tint;
			% peakfl (:, ind) = tpk;   intfl (:, ind) = tint;
		end;
		peakfl (:, ind) = tpk;   intfl (:, ind) = tint;

	end;

	% Print statistics
	fprintf (1, '\nBad fit count: %s\n', num2str(sum (badfit' == 0), ' %2d '));
	fprintf (1, 'Discarded timeslices: %d\n', badtimes);

	% Save the workspace for longer files!
	fsave = strcat (fname, '_save.mat');
	save (fsave, 'callist', 'fit', 'tobs');


	%%%%%%%%%%%%%%%%%%%% Analysis section %%%%%%%%%%%%%%%%%%
	% Plot fitted fluxes
	col = {'b*-', 'm*-', 'r*-', 'k*-', 'g*-', 'y*-', 'w*-', 'c*-'};
	for ind = 1:ncal
		leg_str {ind} = srclist3CR(callist(ind)).name;
	end;

	figure;
	% Peak fluxes and ratios
	subplot (2, 2, 1);
	for ind = 1:ncal
		% pfl = peakfl (ind, badfit(ind, :) == 0);
		% plot (pfl, char(col(ind)));
		plot (peakfl (ind, :), char(col(ind))); % Unfiltered data.
		hold on;
	end;
	title ('Extract peak flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);
	subplot (2, 2, 2);
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
	subplot (2, 2, 3);
	for ind = 1:ncal
		% ifl = intfl (ind, badfit(ind, :) == 0);
		% plot (ifl, char (col(ind)));
		plot (intfl(ind, :), char (col(ind)));
		hold on;
	end;
	title ('Extract int. flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);
	subplot (2, 2, 4);
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
		subplot (2, ncal, ind);
		pfl = peakfl (ind, badfit(ind, :) == 1);
		plotfithist (pfl, 50, gca);
		% hist (pfl, 50);
		title (sprintf ('Peak flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2, ncal, ind+ncal);
		pfl = peakfl (ind, badfitunion == ncal); 
		% hist (pfl./pkflref, 50);
		plotfithist (pfl./pkflref, 50, gca);
		title (sprintf ('Peak flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;

	figure;
	intflref = intfl (1, badfitunion == ncal); 
	for ind = 1:ncal
		subplot (2, ncal, ind);
		ifl = intfl (ind, badfit(ind, :) == 1);
		% hist (ifl, 50);
		plotfithist (ifl, 50, gca);
		title (sprintf ('Integ. flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2,ncal,ind+ncal);
		pfl = intfl (ind, badfitunion == ncal); 
		% hist (pfl./intflref, 50);
		plotfithist (pfl./intflref, 50, gca);
		title (sprintf ('Integ. flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;

	% Generate a plot of the peak position of sources
	xpos = zeros (ncal, ntslices);
	ypos = xpos;
	sigx = xpos; 
	sigy = xpos;
	pa   = xpos;
	for ind=1:ntslices
		for cl = 1:ncal
			fitparam (cl,:) = fit(cl, ind).fitparams;
		end;
		xpos (:,ind) = fitparam (:,1);
		ypos (:,ind) = fitparam (:,2);
		sigx (:,ind) = fitparam (:,4);
		sigy (:,ind) = fitparam (:,5);
		pa (:,ind) = fitparam (:,6);
	end;

	% Plot other fit related data: 
	% histograms of x and y positions, fitted sigx, sigy and position angle
	for ind = 1:ncal
		figure;
		subplot (251);	plot (xpos (ind,:)); title ('lpos (pix)'); axis tight;
		subplot (252);	plot (ypos (ind,:)); title ('mpos (pix)'); axis tight;
		subplot (253);	plot (sigx (ind,:)); title ('sigl (pix)'); axis tight;
		subplot (254);	plot (sigy (ind,:)); title ('sigm (pix)'); axis tight;
		subplot (255);	plot (pa   (ind,:)); title ('pa (rad)'); axis tight;
		
		subplot (256);	plotfithist (xpos (ind,:), 50, gca);
		subplot (257);	plotfithist (ypos (ind,:), 50, gca);
		subplot (258);	plotfithist (sigx (ind,:), 50, gca);
		subplot (259);	plotfithist (sigy (ind,:), 50, gca);
		subplot (2,5,10);	plotfithist (pa   (ind,:), 50, gca);
		mtit (sprintf ('3C%s',srclist3CR(callist(ind)).name), 'yoff',0.025);
	end;
	
