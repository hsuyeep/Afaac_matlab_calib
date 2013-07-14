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
	addpath ../srcfind
	
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
	% NOTE: Should be an all-sky callist which is pruned to the current FoV.
	callist = [200, 133, 237];
	ncal = length (callist);

	% Create datastructures
	tpk = zeros (ncal, 1); tint = zeros (ncal, 1); % Instantaneous peak/int. fl.

	peakfl = zeros (ncal, ntslices);     % Timeseries datastructures
	intfl = peakfl;
	resid = peakfl;

	winsize = 10; % Timeslices           % Window datastructures
	pkwin = zeros (ncal, winsize);
	intwin = pkwin;
	thresh = 5;                         % Rejection threshold = 10*median value

% Debug: Create an image from specified file.
%	img = readimg2bin (fid);
%	map = reshape (img.map, img.pix2laxis, img.pix2maxis);
%   mask = NaN(size (map));
%   mask(meshgrid(img.l).^2 + meshgrid(img.m).'.^2 < 1) = 1;

	% Fill in timeslice window.
	fprintf (1, 'Filling filter window...\n');
	for ind = 1:winsize
		img = readimg2bin (fid);
		[pkwin(:, ind), intwin(:, ind), upcals, resid(:, ind)] = ...
			extractcal (img.map, img.l, img.m, img.tobs, img.freq, rodata, ...
						 srclist3CR, callist, true, 0);	
	end;
	% Generate initial statistics for each calsource, over winsize timeslices.
	% Transpose due to median being a row vector computed over cols.
	pkmed = median (pkwin'); intmed = median (intwin'); 

	% Move back to beginning of file.
	% TODO

	for ind = 1:ntslices
		img = readimg2bin (fid);

		[tpk, tint, upcals, resid(:, ind)] = ...
			extractcal (img.map, img.l, img.m, img.tobs, img.freq, rodata, ...
						 srclist3CR, callist, true, 0);	

		% Decide if current timeslice should enter timeseries.
		% NOTE: Not using the integrated flux!! 
		if (sum (tpk > thresh*pkmed) > 0) % Vectorial comparison.
			fprintf (2, '<-- Discarding time %.2f\n', img.tobs);	
			continue;
		end;

		% Fill window circularly
		winind = mod (ind, winsize) + 1;
		pkwin (:, winind) = tpk; intwin (:, winind) = tint;
		peakfl (:, ind) = tpk; intfl (:, ind) = tint;
	end;

	% Plot fitted fluxes
	col = {'b.-', 'm.-', 'r.-', 'k.-', 'g.-', 'y.-', 'w.-', 'c.-'};
	figure;
	% Peak fluxes and ratios
	subplot (221);
	for ind = 1:ncal
		plot (peakfl(ind, :), char(col(ind)));
		hold on;
	end;
	title ('Extract peak flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	subplot (222);
	for ind = 1:ncal
		plot (peakfl(ind, :)./peakfl(1,:), char(col(ind)));
		% plot (resid(ind,:), char(col(ind)));
		hold on;
	end;
	title ('Extract peak flux ratio of calibrators'); 
	% title ('resid flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');

	% Integrated fluxes and ratios
	subplot (223);
	for ind = 1:ncal
		plot (intfl(ind, :), char (col(ind)));
		hold on;
	end;
	title ('Extract int. flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	subplot (224);
	for ind = 1:ncal
		plot (intfl(ind, :)./intfl(1,:), char(col(ind)));
		hold on;
	end;
	title ('Extract int. flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');

	% Plot histograms
%	figure;
%	for ind = 1:ncal
%		subplot (2,2,ind);
%		hist (peakfl(ind, :));
%	end;
%
%	figure;
%	for ind = 1:ncal
%		subplot (2,2,ind);
%		hist (intfl(ind, :));
%	end;
