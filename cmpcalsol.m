% Script to compare calibration solutions generated on the same data set, but
% using different calibration strategies. Used for characterizing different 
% calibration strategies, or under different observation conditions.
%
% Arguments:
%   nrecs : Number of solution records to compare, -1 for smallest calsol file.
% varargin: list of calib. solution filenames to compare. Expect convcal file
%			name to be the first, as that is used as a reference.	
% pep/07Feb12

function cmpcalsol (nrecs, varargin)

	nfiles = nargin - 1;
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
	col = {'bo', 'mx', 'r+', 'k*', 'gs', 'yd', 'wv', 'c^'};

	% Open all files.
	for ind = 1:nfiles;
		fid(ind) = fopen (varargin{ind}, 'rb'); 
		if (fid(ind) < 0)
			fprintf (2, 'Unable to open file %s\n', varargin{ind});
		end;
	end;

	% Reach the same time instant on all files.
	for ind = 1:nfiles % NOTE: Assuming all files are good!
		sol(ind) = readcalsol (fid(ind));
		tobs(ind) = sol(ind).tobs;
	end;
	t_first = max (tobs);

	for ind = 1:nfiles 
		toff = t_first- sol(ind).tobs + 1;
		for ts = 1:toff
			sol(ind) = readcalsol (fid(ind));
			tobs(ind) = sol(ind).tobs;
	  	end
	end

	% Figure management
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

	% Set up histogram structures: NOTE amp and ph can hold re/im as well!
	ph_err_hist = zeros (length(sol(ind).real_gainsol), histbins);
	% cal_diff = ph_err_hist;
	amp_err_hist = ph_err_hist;
	noise_amp_err_hist = zeros (length(sol(ind).real_sigman), histbins);
	noise_ph_err_hist = noise_amp_err_hist;
	err = zeros (nrecs);

	% Process timeslices
	for ts = 1:nrecs 
		if (sum (isempty(tobs)) ~= 0)
			disp('End of *some* file reached!'); break;
		end;
		toff = sol(ind).tobs - t_first;
		fprintf (1, 'Rec: %04d Time: %.2f, flag ants: %03d.\n', ... 
				 ts, tobs(1), length(sol(1).flagants)); % Assuming same flagged ants

		% Plot LS imaging source fluxes from all calsol files.
		% All sources from a single calsol file are plotted with the same color.
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

		% Percent error between the gain solutions from the two files.
		% NOTE: Assuming there are only two files for now!! FIXME TODO 
		sol1 = complex (sol(1).real_gainsol, sol(1).imag_gainsol);
		sol2 = complex (sol(2).real_gainsol, sol(2).imag_gainsol);
		err (ts) = 100*sum(abs(sol1 - sol2)) / sum(abs(sol1));

		figure (noiseplt);
		plot (err, '-mo');
		xlabel ('Timeslice number');
		ylabel ('Relative percent error');
		title ('Temporal variation of relative error between calib. solutions');

		% Plot gain and phase differences between different solution records.
		figure (gainplt);
		subplot (2,3,1);
		% NOTE: We always use the first calsol file as the reference.
		amp_ts_ref = sol(1).real_gainsol;
		ph_ts_ref = sol(1).imag_gainsol;

		% amp_ts_ref = hypot (sol(1).real_gainsol, sol(1).imag_gainsol);
		% ph_ts_ref = angle (complex (sol(1).real_gainsol, ... 
	    %							sol(1).imag_gainsol));

		noise_amp_ref =  sol(1).real_sigman;
		noise_ph_ref = sol(1).imag_sigman;

		% noise_amp_ref =  hypot (sol(1).real_sigman, sol(1).imag_sigman);
		% noise_ph_ref = angle (complex (sol(ind).real_sigman, ... 
		% 							sol(ind).imag_sigman));

		for ind = 2:nfiles % But currently working on only two...
			amp_ts = sol(ind).real_gainsol;
			% amp_ts = hypot (sol(ind).real_gainsol, sol(ind).imag_gainsol);
			cal_diff = amp_ts_ref - amp_ts; % NOTE: Always assumed 0-mean.
			for ant=1:length(amp_ts) % Generate running histi of all ampdiffs
				bin = int32(histbins/2 + cal_diff(ant)/ampbinwid);
				% Cutoff outliers.
				if (bin > histbins) bin = histbins; end;
				if (bin < 1) bin = 1; end;
				
				amp_err_hist (ant,bin) = amp_err_hist(ant,bin) + 1;
			end;
			plot (ampxaxis, amp_err_hist(43,:), '-o'); % Some random antenna's
			 										   % histogram.
			% hold on;
			plot (ampxaxis, amp_err_hist (80, :), '-ro');
			title ('Calibration solution error histogram');
			xlabel ('Gain solution amplitude error');
			subplot (2,3,2);
			imagesc (amp_err_hist');          % All amp. errors for all antennas.
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
%			title ('Timeseries of gain magnitudes from individual antennas');
		end;	

		subplot (2,3,4);
		for ind = 1:nfiles
			ph_ts = sol(ind).imag_gainsol;
			% ph_ts = angle (complex (sol(ind).real_gainsol, ... 
		  	% 						sol(ind).imag_gainsol));
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
			imagesc (ph_err_hist');          % All amp. errors for all antennas.
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

%{
		figure (noiseplt);
		subplot (2,2,1);
		for ind = 2:nfiles % But currently working on only two...
			noise_amp = hypot (sol(ind).real_sigman, sol(ind).imag_sigman);
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
			noise_ph = angle (complex (sol(ind).real_sigman, ... 
									sol(ind).imag_sigman));
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

		% Read in the next record.
		for ind = 1:nfiles % NOTE: Assuming all files are good!
			sol(ind) = readcalsol (fid(ind));
			tobs(ind) = sol(ind).tobs;
		end;
	end;

	for ind = 1:nfiles
		fclose (fid(ind));
	end;
