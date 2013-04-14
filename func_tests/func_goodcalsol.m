% Script to test the efficiency of judging the goodness of a calibration 
% solution, as determined by comparison with a sliding temporal window of 
% solutions.
% pep/25Mar13
% Arguments:
%  fname : Filename of the calibration solutions obtained via convergent
%			 calibration.
%  offset: Offset in file to begin analysis
%  nrecs : Number of records to operate on.
%  wrbad : bool controlling writing of bad solutions to disk.


function func_goodcalsol (fname, offset, nrecs, wrbad, pltall) 

	debuglev = 4;

	if (debuglev >= 4)
		solparm.gainplt = figure;
		solparm.currsolplt = figure;
		set(0,'Units','pixels') 
		scnsize = get(0,'ScreenSize');
		position = get(solparm.gainplt,'Position');
		outerpos = get(solparm.gainplt,'OuterPosition');
		borders = outerpos - position;
		edge = -borders(1)/2;
		% pos = [left bottom width height]
		pos1 = [edge, 0, scnsize(3)/2 - edge, scnsize(4)/2];
		set(solparm.gainplt,'OuterPosition',pos1);
		pos1 = [edge+scnsize(3)*(1/2), 0, scnsize(3)/2-edge, scnsize(4)/2];
		set(solparm.currsolplt,'OuterPosition',pos1);
	end;

	% Open the calib. solutions file
	fid = fopen (fname, 'rb');

	% Check if we need to write out bad solutions
	if (wrbad == 1)
		k = strfind (fname, '.bin');
		outfname = [fname(1:k-1) '_' num2str(offset) '_badsol.bin'];
		fprintf (1, 'func_goodcalsol:Writing bad sols to file : %s.\n',...
			 outfname);
		fsol = fopen (outfname, 'wb');
	end;

	% Move to offset; not using fseek as records can be of differing sizes!
	for ind=1:offset
		try
			sol0 = readcalsol (fid);
		catch err
			error ('func_goodcalsol: Eof reached!');
		end;
	end;
	Nelem = length (sol0.gainsol);

	if (isempty(nrecs) || nrecs < 0)
		% Determine number of records: Crude way, as recsize could not be 
		% determined correctly!
		t = whos ('sol0');
		recsize = t.bytes; 
		d = dir (fname);
		nrecs = int32 (d.bytes/t.bytes);
		fprintf (1, '-->Filesize: %d, recsize: %d, nrecs: %d\n', ... 
				 d.bytes, t.bytes, nrecs);
	end;

	% Generate data structures 
	% Reject calsolutions with std. across antennas in a station > solthresh
	% percent from mean solutions.
	solparm.solthresh=0.15;
	% solparm.errthresh = 10; % Tolerance as offset from mean gainsol, percent
	% solparm.mse_ph_thresh = 30;% Deg, for mse phase error
	solparm.errthresh = 15; % Tolerance as offset from mean gainsol, percent
	solparm.mse_ph_thresh = 50;% Deg, for mse phase error
	solwinsize = 10; 
	solwindow = complex (zeros (solwinsize, Nelem), ... 
						 zeros (solwinsize, Nelem)); 
	gainmask = zeros (1, Nelem);
	gainmask (sol0.flagant) = 1;
	meangain = zeros (Nelem, nrecs);
	meanerr = meangain;
	stat_std_ts = zeros (6, nrecs);
	badsols = 0;
	col = {'bo-', 'mo-', 'ro-', 'ko-', 'go-', 'yo-', 'wo-', 'co-'};

	% Solution statistics
	histbins = 30;
	re_hist = zeros (Nelem, histbins);
	im_hist = re_hist;
	err_ts = zeros (1, nrecs);
	mse_ph_ts = zeros (1, nrecs);
	mse_amp_ts = mse_ph_ts;
	stnstdhist = zeros (6, histbins); % Histogram of station stds.
	badgainint = 50; % Interval over which to count bad gains, in recs.
	badgainhistory = zeros (1, int32(nrecs/badgainint));
	badint = 1;

	% Fill in the window: 
	fprintf (1, 'Filling solution window...');
	ts = 0;
	while (1)
	% for ts = 1:solwinsize
		try
			sol = readcalsol (fid);
		catch
			fprintf (1, 'func_goodcalsol: Error in readcalsol!\n');
		end;
		% gainsol = complex (sol.real_gainsol, sol.imag_gainsol);
		solstat(1) = goodcalsol (solwindow, ts, sol.gainsol.', gainmask, ...
								 solparm, debuglev);
		if (solstat(1).goodstncal == 1)
			ts = ts + 1;
			if (ts == solwinsize)
				break;
			end;
			solwindow (mod(ts,solwinsize), :) = sol.gainsol;
		end;
	end;
	fprintf (1, 'Done.\n');
	% Mean re/im components for histograms
	meansol0 = mean (solwindow);

	% Generate limits for running histograms
	% Re/im, with 20% ranges.
	rehisthi = max (real(solwindow)) + 0.2*max(real(solwindow));
	rehistlo = min (real(solwindow)) - 0.2*min(real(solwindow));
	imhisthi = max (imag(solwindow)) + 0.2*max(imag(solwindow));
	imhistlo = min (imag(solwindow)) - 0.2*min(imag(solwindow));
	
	% Vector bin widths
%	rebinwid = abs ((rehisthi - rehistlo))/histbins;
%	imbinwid = abs ((imhisthi - imhistlo))/histbins;

	% Scalar bin widths, as some widths were quite close to 0, killing the hists
	rebinwid = max (abs ((rehisthi - rehistlo))/histbins);
	imbinwid = max (abs ((imhisthi - imhistlo))/histbins);
	% xaxis = linspace (rehistlo, rehisthi, histbins);

	% gainsol = complex (sol.real_gainsol', sol.imag_gainsol');
	solstat(1) = goodcalsol (solwindow, solwinsize, sol.gainsol.', ... 
												gainmask, solparm, debuglev);
	% setup histograms for stn. std. devs.
	stnstdhi = solstat(1).stat_std_arr + 1; % max (stnstd);
	stnstdlo = solstat(1).stat_std_arr - 1; %min (stnstd);
	stnbinwid = (stnstdhi - stnstdlo)/histbins;
	% stnxaxis = linspace (stnstdhi, stnstdlo, histbins);

	% Analyze the rest of the timeslices.
	for ts = solwinsize:nrecs
		try
			sol = readcalsol (fid);
		catch
			fprintf (1, 'func_goodcalsol: Error in readcalsol!\n');
		end;
		% gainsol = complex (sol.real_gainsol', sol.imag_gainsol');

		% Operate on the contents of the sliding window.
		meansol = mean(solwindow); % Mean of re/im separately

		% Generate global solution statistics.
		% Gen. running hist of all gainsols, re/im separately
		% NOTE: Can be vectorized with care!
		for ant=1:Nelem 
			% Real part
			bin = int32 (histbins/2 + ...
				(real(sol.gainsol(ant)) - real(meansol0(ant)))./rebinwid);
			% Cutoff outliers.
			if (bin > histbins) bin = histbins; end;
			if (bin < 1) bin = 1; end;
			re_hist (ant,bin) = re_hist (ant,bin) + 1;

			% Imag. part
			bin = int32 (histbins/2 + ...
				(imag (sol.gainsol(ant)) - imag(meansol0(ant)))./imbinwid);
			% Cutoff outliers.
			if (bin > histbins) bin = histbins; end;
			if (bin < 1) bin = 1; end;
			im_hist (ant,bin) = im_hist (ant,bin) + 1;
		end;
		
		% Generate mean re/im. timeseries, mean over sliding window.	
		meangain (:, ts) = meansol;

		% Instantaneous error from mean
		meanerr (:, ts) = meansol - sol.gainsol.';

		% Decide whether to update the solution buffer with current solution.
		solstat(ts) = goodcalsol (solwindow, solwinsize, sol.gainsol.', gainmask, ... 
									solparm, debuglev);
		err_ts(ts) = solstat(ts).err;
		% if (solstat(ts).mse_ph > 2*solparm.mse_ph_thresh)
		%	mse_ph_ts(ts) = solparm.mse_ph_thresh;
		% else
		 	mse_ph_ts(ts) = solstat(ts).mse_ph;
			mse_amp_ts (ts) = solstat(ts).mse_amp;
		% end;
		stat_std_ts (:, ts) = solstat(ts).stat_std_arr;


		% Handle a bad timeslice.
		if (solstat(ts).goodstncal == 0)
			fprintf (2, 'Rec: %03d, Time: %.2f, good = %d.\n', ...
				 ts+offset, sol.tobs, solstat(ts).goodstncal);
			% fprintf (1, '%s\n', num2str(stat_std_arr, '%5.3f '));
			if (wrbad == 1)
				wrcalsol2bin (fsol, sol);
			end;
			badsols = badsols + 1;
		else	
			fprintf (1, 'Rec: %03d, Time: %.2f, good = %d.\n', ...
				 ts+offset, sol.tobs, solstat(ts).goodstncal);
			solwindow (mod(ts,solwinsize) + 1, :) = sol.gainsol;
			% fprintf (1, '%s\n', num2str(stat_std_arr, '%5.3f '));
		end;

		% update std. histograms.
		for stn = 1:6
			bin = int32 (histbins/2 + (solstat(ts).stat_std_arr (stn) - solstat(1).stat_std_arr(stn)) ./ stnbinwid(stn));
			% Cutoff outliers.
			if (bin > histbins) bin = histbins; end;
			if (bin < 1) bin = 1; end;
			stnstdhist (stn, bin) = stnstdhist (stn,bin) + 1;
		end;
		
		if (mod (ts, badgainint) == 0)
			badgainhistory (badint) = badsols;
			badint = badint + 1;
		end;
		% pause;
	end;

	fprintf (1, '\nNumber of bad gains: %d\n', badsols);
	fprintf (1, 'Mean/std of phase std. across stations:\n');
	for stn=1:6
		[m_stn(stn), v_stn(stn), sel]=robustmean (stat_std_ts (stn, :), 3);
		fprintf (1, '%.2f/%.2f  ', m_stn(stn), v_stn(stn));
	end;
	fprintf ('\n');

	% Plot results
	figure;
	subplot (1,2,1);
	% imagesc (re_hist);
	% title ('Histogram of real gainsol');
	hist (err_ts);
	title ('Rel. err. wrt. meansol. histogram');
	[m, v, sel] = robustmean (err_ts, 3);
	fprintf (1, 'Mean/std of rel. err: %f/%f\n', m, v);

	subplot (1,2,2);
	% imagesc (im_hist);
	% title ('Histogram of imag. gainsol');
	hist (mse_ph_ts);
	title (sprintf ('MSE phase histogram, saturating at %d deg.', ... 
					solparm.mse_ph_thresh*2));
	[m, v, sel] = robustmean (mse_ph_ts, 3);
	fprintf (1, 'Mean/std of mse_ph_ts: %f/%f\n', m, v);
	[m, v, sel] = robustmean (mse_amp_ts, 3);
	fprintf (1, 'Mean/std of mse_amp_ts: %f/%f\n', m, v);

	tb1 = uicontrol ('style', 'text');
	set (tb1, 'Units', 'characters');
	pos = get (tb1, 'Position');
    pos(1) = 0; pos (2) = 0; pos(3) = length(fname); pos(4) = 1; 
	set (tb1, 'Position', pos); set (tb1, 'FontSize', 8);
	set (tb1, 'String', fname);
	
	figure;
	subplot (1, 2, 1);
	% imagesc (stnstdhist);
	for stn=1:6
		plot (stnstdhist(stn, :), char(col(stn)));
		hold on;
	end;
	xlabel ('Bin number');
	ylabel ('Counts');
	title ('Histogram of complex std. dev of gain solutions over 6 stations');

	subplot (1, 2, 2);
	plot (badgainhistory, '-ro');
	xlabel (sprintf ('Time in %f rec. units', badgainint));
	ylabel ('Cumulative bad gain count');
	title ('Cumulative bad gain count over time');
	tb2 = uicontrol ('style', 'text');
	set (tb2, 'Units', 'characters');
	pos = get (tb2, 'Position');
    pos(1) = 0; pos (2) = 0; pos(3) = length(fname); pos(4) = 1; 
	set (tb2, 'Position', pos); set (tb2, 'FontSize', 8);
	set (tb2, 'String', fname);

	if (wrbad == 1)
		fclose (fsol);
	end;
	fclose (fid);
	% imagesc (stnstdhist);
