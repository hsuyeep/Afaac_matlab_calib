% Script to compare calibration solutions generated on the same data set, but
% using different calibration strategies. Used for characterizing different 
% calibration strategies, or under different observation conditions.
%
% Arguments:
%   nrecs : Number of solution records to compare, -1 for full length of 
%			smallest calsol file.
%  showplt: Bool, if set to 1, only the relative errors are generated. Faster.
% varargin: list of calib. solution filenames to compare. Expect convcal file
%			name to be the first, as that is used as a reference.	
% pep/07Feb12
%  Added check for timestamp of solution files being compared, to cater to
%  timeslices which can be missed by the tracking calibration.
% pep/17Mar13

function cmpcalsol (nrecs, showplt, varargin)

	nfiles = nargin - 2;
	misstimes = 0;
	histbins = 30;
	ampmean = 0; phmean = 0;
	amphisthi = 0.5; amphistlo = -0.5;
	ampbinwid = (amphisthi - amphistlo)/histbins;
	ampxaxis = linspace (amphistlo, amphisthi, histbins);
	phhisthi =  0.5; phhistlo = -0.5; % Rad, if phase.
	phbinwid = (phhisthi - phhistlo)/histbins;
	phxaxis = linspace (phhistlo, phhisthi, histbins);
	fprintf (1, 'Total number of files specified: %d.\n', nfiles);
	fid  = zeros (1, nfiles);

	% Indicates waiting a timeslice as a nop, for temporal alignment of sols. 
	% If waitind (ind) == 1 for a given file, at a given time, ignore it.
	waitind = zeros (1, nfiles); 
	col = {'bo', 'mx', 'r+', 'k*', 'gs', 'yd', 'wv', 'c^'};

	% Open all files.
	for ind = 1:nfiles;
		fprintf (1, 'Opening file: %s at index %d\n', varargin{ind}, ind);
		fid(ind) = fopen (varargin{ind}, 'rb'); 
		if (fid(ind) < 0)
			fprintf (2, 'Unable to open file %s\n', varargin{ind});
			return;
		end;
	end;

	% Determine number of records to operate on
	if (nrecs == -1)
		for ind = 1:nfiles
			% Determine number of records: Crude way, as recsize could not be 
			% determined correctly!
			try
				rec0 = readcalsol (fid(ind));
			catch err
				fprintf (2, 'cmpcalsol: Eof reached!\n');
				return;
			end;
			t = whos ('rec0');
			recsize = t.bytes; 
			d = dir (varargin{1});
			ntimes = int32 (d.bytes/t.bytes);
			fprintf (1, '-->Filesize: %d, recsize: %d, nrecs: %d\n', ... 
				 d.bytes, t.bytes, ntimes);
			fsize (ind) = ntimes;
		end;
	else
		fsize = nrecs*ones (1,nfiles);
	end;
	nrecs = int32(min (fsize));
	fprintf (1, 'Setting number of records to %d\n', nrecs);

	% Reach the same time instant on all files.
	for ind = 1:nfiles % NOTE: Assuming all files are good!
		sol(ind) = readcalsol (fid(ind));
		tobs(ind) = sol(ind).tobs;
	end;
	t_first = max (tobs);
	sol_first = sol(1); % NOTE: Gain deviations wrt. first time record!.

	for ind = 1:nfiles 
		toff = t_first- sol(ind).tobs + 1;
		for ts = 1:toff
			sol(ind) = readcalsol (fid(ind));
			tobs(ind) = sol(ind).tobs;
	  	end
	end

	% Figure management
	if (showplt == 1)
    	noiseplt = figure;
    	set (gca, 'Fontsize', 14);
	else
    	fluxplt = figure;
    	set (gca, 'Fontsize', 14);
    	gainplt = figure;
    	set (gca, 'Fontsize', 14);
    	noiseplt = figure;
    	set (gca, 'Fontsize', 14);
    	set(0,'Units','pixels') 
    	scnsize = get(0,'ScreenSize');
    
    	position = get(fluxplt,'Position');
    	outerpos = get(fluxplt,'OuterPosition');
    	borders = outerpos - position;
    	edge = -borders(1)/2;
    	% pos = [left bottom width height]
    	pos1 = [edge, scnsize(4) * (1/2), scnsize(3)/2 - edge, scnsize(4)/2];
    	set(fluxplt,'OuterPosition',pos1);
    	pos1 = [edge, 0, scnsize(3) - edge, scnsize(4)/2];
    	set(gainplt,'OuterPosition',pos1);
    	pos1 = [edge+scnsize(3)/2, scnsize(4) * (1/2), scnsize(3)/2 - edge, ... 
    			scnsize(4)/2];
    	set(noiseplt,'OuterPosition',pos1);
	end;

	% Set up histogram structures: NOTE amp and ph can hold re/im as well!
	ph_err_hist = zeros (length(sol(ind).gainsol), histbins);
	% cal_diff = ph_err_hist;
	amp_err_hist = ph_err_hist;
	noise_amp_err_hist = zeros (length(sol(ind).sigman), histbins);
	noise_ph_err_hist = noise_amp_err_hist;
	err = zeros (1, nrecs);

	% Generate log file for plotting results
	fname = sprintf ('cmpcalsol_%s.txt', datestr (now));
	flog = fopen (fname, 'wt');
	fprintf (flog, '%% Gain solution difference between files %s and %s.\n',...
			 varargin{1}, varargin{2});		    
	fprintf (flog, '%% Generated by cmpcalsol.m at %s.\n', datestr (now));
	fprintf (flog, '%% Time (MJDSec) fluxdiff phaseerr\n');
	
	% Process timeslices
	for ts = 1:nrecs 
		if (sum (isempty(tobs)) ~= 0)
			disp('End of *some* file reached!'); break;
		end;
		toff = sol(ind).tobs - t_first;
		fprintf (1, '\nRec: %04d, Time:', ts);
		for ind=1:nfiles 
			fprintf (1, '%.2f/%d ', tobs(ind), length(sol(ind).flagant));
		end;

		% Plot LS imaging source fluxes from all calsol files.
		% All sources from a single calsol file are plotted with the same color.
	 	if (sum (waitind == 0) >= 2) % At least 2 files needed for comparison!
			if (showplt == 0)
				figure (fluxplt);
				set (gca, 'Fontsize', 14);
				for ind = 2:nfiles
					for src = 1:sol(ind).calsrcs
						plot (toff, sol(1).sigmas(src) - sol(ind).sigmas(src), ... 
							  char(col(src)));
							 %  char(col(ind)));
						hold on;
					end;
				end;
				title (sprintf ('Extracted fluxes of model sources.'));
				xlabel ('Time (sec from obs. commencement)');
				ylabel ('Flux ratio normalized to CasA flux');
			end;

			% Percent error between the gain solutions from the two files.
			% NOTE: Assuming there are only two files for now!! FIXME TODO 
			flag1 = zeros (1, 288); flag1 (sol(1).flagant) = 1;
			srcsel = (sol(1).sigmas ~= 0);
			% flag2 = zeros (1, 288); flag2 (sol(2).flagant) = 1;
			
			% Solutions with complex gains
%			sol1 = [sol(1).gainsol(flag1 == 0); sol(1).sigmas(srcsel); sol(1).thsrc_wsf(srcsel); sol(1).phisrc_wsf(srcsel)];
%			sol2 = [sol(2).gainsol(flag1 == 0); sol(2).sigmas(srcsel); sol(2).thsrc_wsf(srcsel); sol(2).phisrc_wsf(srcsel)];
%			% err (ts) = 100*sum(abs(sol1 - sol2)) / sum(abs(sol1));
%			err (ts) = (100/length(sol1))*sum( abs(sol1 - sol2) ./ abs(sol1));

			% Solutions with explicit real and imaginary gains.
%			sol1 = [real(sol(1).gainsol(flag1 == 0)); imag(sol(1).gainsol(flag1 == 0)); sol(1).sigmas(srcsel); sol(1).thsrc_wsf(srcsel); sol(1).phisrc_wsf(srcsel)];
%			sol2 = [real(sol(2).gainsol(flag1 == 0)); imag(sol(2).gainsol(flag1 == 0)); sol(2).sigmas(srcsel); sol(2).thsrc_wsf(srcsel); sol(2).phisrc_wsf(srcsel)];

			% err (ts) = 100*sum(abs(sol1 - sol2)) / sum(abs(sol1));
			% err (ts) = (100/length(sol1))*sum( abs(sol1 - sol2) ./ abs(sol1));
			sol1 = angle (sol(1).gainsol(flag1 == 0)); 
			sol2 = angle (sol(2).gainsol(flag1 == 0)); 
			phdiff = sol1-sol2;
			err (ts) = mean( abs(phdiff(phdiff~=0))); % Mean Abs. Deviation.
			% err (ts) = sqrt (mean((phdiff(phdiff~=0).^2))); % Mean Abs. Deviation.
			fprintf (flog, '%f %f\n', sol(ind).tobs, err(ts));

			% err (ts) = 100*sum(abs(sol1 - sol2)) / sum(abs(sol1));
			fprintf (1, 'sol1:%f, sol2:%f, Err: %f.\n', sol1, sol2, err(ts));
	
			figure (noiseplt);
			plot (toff, err(ts), '-ro');
			hold on;
			xlabel ('Timeslice number');
			ylabel ('Relative percent error');
			title ('Temporal variation of relative error between calib. sols.');
	
			if (showplt == 0)
				% Plot gain and phase differences between different sol. records.
				figure (gainplt);
				subplot (2,3,1);
				% NOTE: We always use the first calsol file as the reference.
				amp_ts_ref = real (sol(1).gainsol);
				ph_ts_ref = imag (sol(1).gainsol);
		
				% amp_ts_ref = abs (sol(1).gainsol, sol(1).gainsol);
				% ph_ts_ref = angle (sol(1).gainsol)
		
				noise_amp_ref =  real (sol(1).sigman(:));
				noise_ph_ref = imag (sol(1).sigman(:));
		
				% noise_amp_ref =  hypot (sol(1).sigman);
				% noise_ph_ref = angle (sol(1).sigman);
		
				for ind = 2:nfiles % But currently working on only two...
					amp_ts = real (sol(ind).gainsol);
					% amp_ts = abs (sol(ind).gainsol);
					cal_diff = amp_ts_ref - amp_ts; % NOTE: Always assumed 0-mean.
					for ant=1:length(amp_ts) % Gen. running histi of all ampdiffs
						bin = int32(histbins/2 + cal_diff(ant)/ampbinwid);
						% Cutoff outliers.
						if (bin > histbins) bin = histbins; end;
						if (bin < 1) bin = 1; end;
						
						amp_err_hist (ant,bin) = amp_err_hist(ant,bin) + 1;
					end;
					plot (ampxaxis, amp_err_hist(43,:), '-o'); % Some random ant's
					 										   % histogram.
					% hold on;
					plot (ampxaxis, amp_err_hist (80, :), '-ro');
					title ('Calibration solution error histogram');
					xlabel ('Gain solution amplitude error');
					subplot (2,3,2);
					imagesc (amp_err_hist');     % All amp. errors for all antennas.
					title ('Calibration solution error histogram');
					xlabel ('Gain solution amplitude error');
					
					subplot (2,3,3);
					plot (cal_diff, 'o');
					xlabel ('Antenna number');
					ylabel ('Gain amplitude difference');
		
		%			for stat = 1:6 % 6 stations
		%				plot (toff, amp_ts_ref(24*stat) - amp_ts(24*stat), ... 
		%					  char(col(ind)));
		%				hold on;
		%			end;
		%			title ('Timeseries of gain magnitudes from individual ants');
				end;	
		
				subplot (2,3,4);
				for ind = 1:nfiles
					ph_ts = imag (sol(ind).gainsol);
					% ph_ts = angle (sol(ind).gainsol);
					cal_diff = ph_ts_ref - ph_ts; % NOTE: Always assumed 0-mean.
					for ant=1:length(ph_ts)
						bin = int32(histbins/2 + cal_diff(ant)/phbinwid);
						% Cutoff outliers.
						if (bin > histbins) bin = histbins; end;
						if (bin < 1) bin = 1; end;
						ph_err_hist (ant,bin) = ph_err_hist(ant,bin) + 1;
					end;
					plot (phxaxis, ph_err_hist (43, :), '-o');
					% hold on;
					plot (phxaxis, ph_err_hist (80, :), '-ro');
					xlabel ('Gain solution phase error');
					subplot (2,3,5);
					imagesc (ph_err_hist');      % All amp. errors for all antennas.
					xlabel ('Gain solution phase error');
		
					subplot (2,3,6);
					plot (cal_diff, 'o');
					xlabel ('Antenna number');
					ylabel ('Gain phase difference (rad)');
		%			for stat = 1:6 % 6 stations
		%				plot (toff, ph_ts_ref (24*stat) - ph_ts(24*stat), ... 
		%					  char(col(ind)));
		%				hold on;
		%			end;
		%			title ('Timeseries of gain phases from individual antennas');
				end;
			end;
		end;

%{
		figure (noiseplt);
		subplot (2,2,1);
		for ind = 2:nfiles % But currently working on only two...
			noise_amp = abs (sol(ind).sigman);
			noise_diff = noise_amp_ref - noise_amp; % NOTE: Always assumed 0-mean.
			for ant=1:length(noise_amp) % Generate running histi of all ampdiffs
				bin = int32(histbins/2 + noise_diff(ant)/ampbinwid);
				% Cutoff outliers.
				if (bin > histbins) bin = histbins; end;
				if (bin < 1) bin = 1; end;
				
				noise_amp_err_hist (ant,bin) = noise_amp_err_hist(ant,bin) + 1;
			end;
			plot (ampxaxis, noise_amp_err_hist(43,:), '-o');
			title ('Calibration solution error histogram');
			xlabel ('Gain solution amplitude error');
			
			subplot (2,2,2);
			plot (noise_diff, 'o');
			xlabel ('Antenna number');
			ylabel ('Gain amplitude difference');

%			for stat = 1:6 % 6 stations
%				plot (toff, noise_amp_ref(24*stat) - noise_amp(24*stat), ... 
%					  char(col(ind)));
%				hold on;
%			end;
%			title ('Timeseries of gain magnitudes from individual antennas');
		end;	

		subplot (2,2,3);
		for ind = 1:nfiles
			noise_ph = angle (sol(ind).sigman);
			noise_diff = noise_ph_ref - noise_ph; % NOTE: Always assumed 0-mean.
			for ant=1:length(noise_ph)
				bin = int32(histbins/2 + noise_diff(ant)/phbinwid);
				% Cutoff outliers.
				if (bin > histbins) bin = histbins; end;
				if (bin < 1) bin = 1; end;
				noise_ph_err_hist (ant,bin) = noise_ph_err_hist(ant,bin) + 1;
			end;
			plot (phxaxis, noise_ph_err_hist (43, :), '-o');
			xlabel ('Gain solution phase error');

			subplot (2,2,4);
			plot (noise_diff, 'o');
			xlabel ('Antenna number');
			ylabel ('Gain phase difference');
%			for stat = 1:6 % 6 stations
%				plot (toff, noise_ph_ref (24*stat) - noise_ph(24*stat), ... 
%					  char(col(ind)));
%				hold on;
%			end;
%			title ('Timeseries of gain phases from individual antennas');
		end;
%}

		% Read in the next record of the reference file..
		try
			sol(1) = readcalsol (fid(1));
		catch
			fprintf (2, 'cmpcalsol: Error in reading file, quitting.\n');
			fprintf (1, 'Total missing times between files: %d.\n', misstimes);

			for nfile = 1:nfiles
				fclose (fid(nfile));
			end;
		end;

		tobs(1) = sol(1).tobs;
		for ind = 2:nfiles % NOTE: Assuming all files are good to read!
			% Only read if we're not waiting.
			if (waitind(ind) == 0) 
				try
					sol(ind) = readcalsol (fid(ind));
				catch
					fprintf (2,'cmpcalsol: Error in reading file, quitting.\n');
					fprintf (1,'Total missing times between files: %d.\n', ...
							 misstimes);
					for nfile = 1:nfiles
						fclose (fid(nfile));
						break;
					end;
				end;
			end;
			tobs(ind) = sol(ind).tobs;
			waitind (ind) = 0; % Will be reset to 1 if there is still an offset.

			if (tobs(ind) ~= tobs(1)) % NOTE: Always wrt. first file!
				misstimes = misstimes + 1;
				offset = int32 (tobs(1) - tobs(ind));
				fprintf (2, '\n<--Missing Time:%.2f   %.2f, offset: %f.', ...
						 tobs(1), tobs(ind), offset);
				if (offset < 0)
					% Moved ahead of ref. file, need to wait.
					waitind (ind) = 1;
				else
					% Move till we reach the ref. position.
					fprintf (2, 'Skipping %d recs on file %d.', offset, ind);
					for ts=1:offset
						try
							sol(ind) = readcalsol (fid(ind));
						catch 
							fprintf (2, 'cmpcalsol: EoF on sol file!\n');
							break;
						end;
						% NOTE: We may not have to move all 'offset' records, as
						% we may be missing records ourselves!
						if (sol(ind).tobs >= tobs(1))
							break;
						end;
					end;
					tobs(ind) = sol(ind).tobs;
				end;
			end;
		end;
		% pause; 
	end;
	fprintf (1, 'Total missing times between files: %d.\n', misstimes);

	for ind = 1:nfiles
		if (fid (ind) > 0)
			fclose (fid(ind));
		end;
	end;

   if (flog) fclose (flog); end;
