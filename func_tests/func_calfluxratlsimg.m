% Script to excercise xcallsimg.m
% pep/17Jul13

function [fit, gainsolsum] = func_calfluxratlsimg (fname, fcalsol, offset, ntslices, posfilename, debug)

	% Calibrated visibilities
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('func_calfluxrat: fid < 0! Quitting.');
		return;
	end;
		
	% Load calibration solution file for solution timeseries.
	fsol = fopen (fcalsol, 'rb');
	try
		rec0 = readcalsol (fsol);
	catch err
		fprintf (2, 'func_calfluxratanal: Eof reached!\n');
		return;
	end;
	rec = readcalsol (fsol);
	dt = rec.tobs - rec0.tobs;
	flagant = zeros (1, 288);
	flagant (rec.flagant) = 1;
	
	fseek (fid, 0, 'bof');

	% Load various meta data
    rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
    rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
    rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	load (posfilename, 'posITRF', 'poslocal');
    rodata.posITRF = posITRF; 	% NOTE: Not removing flagged vis. for now.
    rodata.poslocal = poslocal;
	load srclist3CR.mat

	if (ntslices < 0)
		[ntslices, tmin, tmax, dt] = getnrecs (fname);
	end;

	% NOTE: Current list assumes 3CR catalog!
	% calsrcs = 3C[295, 219, 338, 410.1, 10, 84, 338]
	% callist = [200, 117, 133, 237];  % For day time observations
	callist = [4, 54, 314, 265, 256];% For nite time observations
	ncal = length (callist);
	rodata.srcsel = callist;         % For compatibility with wsf_srcpos.m

	fit     = zeros (ncal,  ntslices);
	fitth   = zeros (ncal,  ntslices);
	fitphi  = zeros (ncal,  ntslices);
	tobs    = zeros (ntslices);
	re_gain = zeros (ntslices, 288);
	im_gain = zeros (ntslices, 288);
	acorr   = zeros (ntslices, 288);
	acorrts = zeros (ntslices, 1);
	srcflux = zeros (ntslices, rec.calsrcs);
	sysnoise= zeros (ntslices, 288);

	for ind = 1:ntslices
		% Read in calibrated visibilities
		try
			[acc, tobs(ind), freq] = readms2float (fid, -1, -1, 288);
		catch err
			fprintf (2, 'func_calfluxratlsimg: Eof of calvis file reached!\n');
			break;
		end;

		% Read in cal. solutions
		rec = readcalsol (fsol);
		if (rec.tobs ~= tobs (ind))
			fprintf (1, '###  Calsol-img time offset: %f at ind %d!\n', ...
					rec.tobs-tobs(ind), ind); 
		end;

		% Read in a corresponding calsol
		if (rec.tobs < tobs(ind))
			while (rec.tobs < tobs(ind))
				ind = ind + 1;
				rec = readcalsol (fsol);
				re_gain (ind, :) = real (rec.gainsol);
				im_gain (ind, :) = imag (rec.gainsol);
			end;
		else
			while (tobs(ind) < rec.tobs)
				ind = ind + 1;
				[acc, tobs(ind), freq] = readms2float (fid, -1, -1, 288);
			end;
		end;
		
		flagant = zeros (1, 288);
		flagant (rec.flagant) = 1;
		% ---- In case of flagged antennas, generate reshaped ACM ---- 
	    if (length(rec.flagant) ~= 0)
	        disp (['Flagant: ' num2str(rec.flagant')]);
	        antmask = zeros (size (acc));
	        posmask = zeros (size (rodata.posITRF));
	        rem_ants = length(acc) - length(rec.flagant);
	        for fl = 1:length(rec.flagant)
	            antmask (rec.flagant(fl), :) = 1; 
				antmask (:,rec.flagant(fl)) = 1;
	            posmask (rec.flagant(fl), :) = 1;
	        end
	        acc = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
	        rodata.posITRF_fl = reshape(rodata.posITRF(posmask ~=1), ...
										[rem_ants, 3]);
	        disp (['NOTE: ACM resized after flagging to ' num2str(size (acc))]);
	    else
	        rodata.posITRF_fl = rodata.posITRF;
	    end

		re_gain(ind, :) = real (rec.gainsol);
		im_gain(ind, :) = imag (rec.gainsol);
		srcflux (ind, :) = rec.sigmas';
		sysnoise (ind, flagant == 0) = diag (rec.sigman);
		fprintf (1, '\nRec %4d ', ind);
		[upcals, fittemp] = ...
			xcallsimg (acc, tobs(ind), freq, rodata, ...
						 srclist3CR, callist, true, 0);	
		fit (:,ind) = fittemp.peakfl;	
		fitth (:, ind) = fittemp.thsrc_wsf;
		fitphi(:, ind) = fittemp.phisrc_wsf;
		% fprintf (1, '%s\n', num2str(fit.peakfl, '%f '));
	end;
	
	% Plot fluxes
	col = {'b*-', 'm*-', 'r*-', 'k*-', 'g*-', 'y*-', 'w*-', 'c*-'};
	for ind = 1:ncal
		leg_str {ind} = srclist3CR(callist(ind)).name;
	end;

	gainsol = hypot (re_gain, im_gain);
	gainsolsum = sum (gainsol');
	figure;
	subplot (611);
	plot (gainsolsum); grid on; axis tight;
	title ('Average gain solution amp. timeseries');
	for ind = 1:ncal
		subplot (6, 1, ind+1);
		plot (conv (fit (ind, :), ones (1,10), 'same'),  char(col(ind)));
		grid on; axis tight;
		title (sprintf ('3C%s', srclist3CR(callist(ind)).name));
	end;
	title ('Extract peak flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	% legend (leg_str);

%	for ind = 1:ncal
%		plot (fit(ind, :)./fit(1,:), char(col(ind)));
%		% plot (resid(ind,:), char(col(ind)));
%		hold on;
%	end;
%	title ('Extract peak flux ratio of calibrators'); 
	% title ('resid flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	%legend (leg_str);

	% Plot histograms. 
	figure;
	for ind = 1:ncal
		subplot (2,5,ind);
		hist (fit(ind,:));
		title (sprintf ('Peak flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2,5,ind+ncal);
		hist (fit(ind,:)./fit(1,:));
		title (sprintf ('Peak flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;

	% Compute correlation coefficients between sources
	fprintf (1, 'Corr.coe between monitor sources and ave. gainsol (last col).\n');
	[r, p] = corrcoef ([fit' gainsolsum']) 
	
	fclose (fid);
	fclose (fsol);
