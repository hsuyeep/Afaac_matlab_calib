% Script to generate a plot of the timeseries of the norm of visibilities 
% provided as input, and its 1st derivative.
% pep/18Oct12
% Updated to produce a histogram of norm values, in order to detect outliers.
% pep/29Jan13

function [times, normts] = gentsnorm (fname, ntimes, showplt)
	fid = fopen (fname, 'rb');
 	if ntimes == -1
		[ntimes, tmin, tmax, dt] = getnrecs (fname);
	end;

	winsize = 10;
	window = zeros (1, winsize);

	normts = zeros (1, ntimes);
	movavg = zeros (1, ntimes - winsize - 1); % Ignore first window size values.
	movmed = zeros (1, ntimes - winsize - 1);
	movmad = zeros (1, ntimes - winsize - 1);
	medianthresh = 2.5; % Reject timeslices with median > medianthresh*movmed;

	normtsdiff = zeros (1, ntimes);
	times = zeros (1, ntimes);

	% Fill the window
	for ind = 1:winsize
		[acc, times(ind), freq] = readms2float (fid, -1, -1, 288);
		acc (isnan(acc)) = 0; % NOTE NOTE! Norm fails if members are NaNs!
		window(ind) = norm (acc, 'fro');
		normts(ind) = window(ind);
	end;

	% Find the mean/median/mad over the window
	movavg (ind) = mean (window);
	movmed (ind) = median (window);
	movmad (ind) = median (abs (window - movmed (ind)));
	

	statind = 1; % Separate counter for statistics.
	for ts=ind+1:ntimes
		[acc, times(ts), freq] = readms2float (fid, -1, -1, 288);

%		% Skipping bad data stretches, for LBA_OUTER_BAND60/*4min*cal.bin
%		if (ts == 35 || ts == 68 || ts == 80 || ts == 151 || ts == 234)
%			disp(sprintf ('Skipping timeslice %d, time %f.', ts, times(ts)));
%			normts(ts) = normts(ts-1); % Just to make it not go to 0;
%			normtsdiff(ts-1) = normtsdiff(ts-2); % Just to make it not go to 0;
%			continue;
%		end
		acc (isnan(acc)) = 0; % NOTE NOTE! Norm fails if members are NaNs!
		normts(ts) = norm (acc, 'fro');
		normtsdiff(ts-1) = normts(ts) - normts(ts-1);

		window (mod (ts, winsize) + 1) = normts (ts);
		movavg (statind) = mean (window);
		movmed (statind) = median (window);
		movmad (statind) = median (abs (window - movmed (ind)));

		% Check for bad timeslices.
		if (normts(ts) > medianthresh*movmed(statind))
			fprintf (2, '<-- Discarding rec: %03d, Time: %.2f. Excess: %.2f\n', ...
				  ts, times (ts), normts(ts)./(medianthresh*movmed(statind)));
			continue;
		end;
		statind = statind + 1;
		% if (mod (ts, 100) == 99) disp ('.');
	end;

	if (showplt ~= 0)
	    t_first = times(1);
		num = mjdsec2datenum (t_first);
		subplot (2, 3, 1);
		% plot (times, normts);
		plot ((times-t_first)/86400+num, normts);
		grid on; axis tight;
		title (sprintf ('Frobenius norm time series.'));
		datetick ('x', 13, 'keepticks'); % Print HH:MM:SS legend on the time axis.	
	
		subplot (2, 3, 2);
		plot ((times-t_first)/86400+num, normtsdiff);
		grid on; axis tight;
		title (sprintf ('Frobenius norm time series 1st derivative.'));
		datetick ('x', 13, 'keepticks'); % Print HH:MM:SS legend on the time axis.	
	
		subplot (2, 3, 3);
		hist (normts, 100);
		title (sprintf ('Frobenius norm histogram.'));
	
		subplot (2, 3, 4);
		plot ((times(winsize+1:end)-t_first)/86400+num, movavg);
		grid on; axis tight;
		title (sprintf ('Frobenius norm moving average.'));
		datetick ('x', 13, 'keepticks'); % Print HH:MM:SS legend on the time axis.	
	
		subplot (2, 3, 5);
		plot ((times(winsize+1:end)-t_first)/86400+num, movmed);
		grid on; axis tight;
		title (sprintf ('Frobenius norm moving median.'));
		datetick ('x', 13, 'keepticks'); % Print HH:MM:SS legend on the time axis.	
	
		subplot (2, 3, 6);
		plot ((times(winsize+1:end)-t_first)/86400+num,movmad);
		grid on; axis tight;
		title (sprintf ('Frobenius norm moving MAD.'));
		datetick ('x', 13, 'keepticks'); % Print HH:MM:SS legend on the time axis.	
	end;
