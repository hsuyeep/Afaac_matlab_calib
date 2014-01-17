% Script to carry out analysis of source fluxes extracted using func_calfluxrat
% pep/23Jul13
% Arguments:
%	winsize: Size of smoothening tophat kernel to apply.
%   recoff : record offset from which to start processing
%  ntslices: Number of consecutive timeslices to consider. -1 indicates all.
%   order  : Order of polynomial to fit on flux timeseries for each source.

function func_calfluxratanal (winsize, recoff, ntslices, order)
	% Night time data
	% Load extracted flux timeseries (generated after running func_calfluxrat.m)
	load ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB002_LBA_OUTER_SPREAD_1ch_1_convcal_1_el_fftimg.bin_save.mat');
	% Load calibration solution file for solution timeseries.
	fid = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB002_LBA_OUTER_SPREAD_1ch_1_convcalsol.bin', 'rb');
	% Load uncalibrated visibility file for pulling out autocorrelations.
	funcal = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch.bin', 'rb');
	load '../srclist3CR.mat';

	% Day time data.
	% Load extracted flux timeseries
%	load ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_1_convcal_1_el_fftimg.bin_save.mat');
%	% Load calibration solution file for solution timeseries.
%	fid = fopen ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_1_convcalsol.bin', 'rb');
%	% Load uncalibrated visibility file for pulling out autocorrelations.
%	funcal = fopen ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr.bin', 'rb');
	

	try
		rec0 = readcalsol (fid);
	catch err
		fprintf (2, 'pltcalsol: Eof reached!\n');
		return;
	end;
	rec = readcalsol (fid);
	dt = rec.tobs - rec0.tobs;
	flagant = zeros (1, 288);
	flagant (rec.flagant) = 1;
	
	fseek (fid, 0, 'bof');
	% Uncalibrated files always lag due to the windowing during calibration
	[acc, tacc, freq] = readms2float (funcal, -1, -1, 288);
	while (tacc < rec0.tobs)
		[acc, tacc, freq] = readms2float (funcal, -1, -1, 288);
	end;
	fseek (funcal, -8*(2+(288*(288+1)/2)), 'cof'); % Move back one record.

	% Figure out number of sources to be monitored, and over what timescale.
	[ncal totntslices] = size (fit);
	if (ntslices < 0) ntslices = totntslices - recoff; end;
	fprintf (1, 'Sources: %d, Timeslices (Total/processed): %d/%d\n', ...
			 ncal, totntslices, ntslices);
	recend = recoff + ntslices;

	% Figure management
	numxticks = 4;
	peakplt = figure;
	posplt = figure;
%	set(0,'Units','pixels') 
%	scnsize = get(0,'ScreenSize');
%	
%	position = get(posplt,'Position');
%	outerpos = get(posplt,'OuterPosition');
%	borders = outerpos - position;
%	edge = -borders(1)/2;
%		% pos = [left bottom width height]
%	pos1 = [edge, 0, scnsize(3) - edge, scnsize(4)];
%	set(posplt,'OuterPosition',pos1);

	% Generate required timeseries
	peakfl = zeros (ncal, ntslices);
	intfl = peakfl;
	xpos = peakfl;
	ypos = xpos;
	sigx = xpos; 
	sigy = xpos;
	% tobs = zeros (ntslices, 1);
	re_gain = zeros (ntslices, 288);
	im_gain = zeros (ntslices, 288);
	acorr   = zeros (ntslices, 288);
	acorrts = zeros (ntslices, 1);
	srcflux = zeros (ntslices, rec.calsrcs);
	sysnoise= zeros (ntslices, 288);

	ts = 1; % local index for dealing with record offsets.
	for ind = recoff:recend - 1

		%%%% Handle calibration solution 
		rec = readcalsol (fid);
		if (rec.tobs ~= tobs (ind))
			fprintf (1, '###  Calsol-img time offset: %f at ind %d!\n', ...
					rec.tobs-tobs(ind), ind); 
		end;
		% Read in a corresponding calsol
		while (rec.tobs < tobs(ind))
			rec = readcalsol (fid);
		end;
		flagant = zeros (1, 288);
		flagant (rec.flagant) = 1;

		re_gain(ts, :) = real (rec.gainsol);
		im_gain(ts, :) = imag (rec.gainsol);
		srcflux (ts, :) = rec.sigmas';
		sysnoise (ts, flagant == 0) = diag (rec.sigman);

		%%% Handle autocorrelation.
		[acc, acorrts(ts), freq] = readms2float (funcal, -1, -1, 288);
		if (acorrts(ts) ~= tobs (ind))
			fprintf (1, '###  Autocorr-img time offset: %f at ind %d!\n', ...
					acorrts(ts)-tobs(ind), ind); 
			while (tacc < rec.tobs)
				[acc, tacc, freq] = readms2float (funcal, -1, -1, 288);
			end;
		end;
		acorr (ts, :) = diag (acc)';

		%%% Handle monitored sources.
		for cl = 1:ncal
			fitparam (cl,:) = fit(cl, ind).fitparams;
		end;
		tpk = fitparam (:,3);
		tint = 2*pi*fitparam (:,3).*fitparam (:,4).*fitparam (:,5); 
		peakfl(:, ts) = tpk;
		intfl(:, ts) = tint;
		xpos (:,ts) = fitparam (:,1);
		ypos (:,ts) = fitparam (:,2);
		sigx (:,ts) = fitparam (:,4);
		sigy (:,ts) = fitparam (:,5);
		fprintf (1, 'ind: %d flagant: %d, sigman size : %d %d, diag: %d\n', ind, sum(rec.flagant == 1), size (rec.sigman, 1), size(rec.sigman, 2), length (diag(rec.sigman)));

		ts = ts + 1; % increment local time counter
	end;

	col = {'b*', 'm*', 'r*', 'k*', 'g*', 'y*', 'w*', 'c*'};

	% Show autocorr timeseries
	figure;
	for ind = 1:6 % Plot for 6 stations
		gain = conv (acorr (:, 24*ind), 1/winsize*ones (1, winsize), 'same');
		plot (gain,  char(col(ind)));
		hold on;
	end;
	title ('Uncalibrated autocorrelation timeseries from 6 antennas (24*ind)');
	xlabel ('Timeslice number'); ylabel ('Autocorrelation');
	grid;
	% Show Correlation between autocorrelations
	[r_acorr, p_acorr] = corrcoef (acorr);
	figure; imagesc (r_acorr); colorbar;
	xlabel ('Ant. number'); ylabel ('Antenna number');
	title ('Corr.coe between Autocorrs of all antennas');

	figure;
	for ind = 1:6 % Plot for 6 stations
		gain = conv (abs(sysnoise (:, 24*ind)), 1/winsize*ones (1, winsize), 'same');
		plot (gain,  char(col(ind)));
		hold on;
	end;
	title ('Estimated system noise timeseries from 6 antennas (24*ind)');
	xlabel ('Timeslice number'); ylabel ('System noise (arb. units)');
	grid;
	% Show Correlation between noise timeseries
	[r_acorr, p_acorr] = corrcoef (abs(sysnoise));
	figure; imagesc (r_acorr); colorbar;
	xlabel ('Ant. number'); ylabel ('Antenna number');
	title ('Corr.coe between system noise of all antennas');

	figure;
	% Show model source flux timeseries
	for ind = 1:rec.calsrcs
		plot (conv (srcflux(:, ind), 1/winsize*ones (1, winsize), 'same'), char (col(ind)));
		hold on;
	end;
	title ('Model source LS imaging flux estimate timeseries');
	xlabel ('Timeslices'); ylabel ('Flux (normalized to CasA)');
	legend ('3C461 (CasA)', '3C405 (CygA)', '3C144 (TauA)', '3C274 (VirA)', 'Sun');
	grid;

	% Show calib gain sols.
	figure;
	for ind = 1:6 % Plot for 6 stations
		gain = conv (abs(complex (re_gain (:, 24*ind), im_gain (:, 24*ind))),...
					 1/winsize*ones (1, winsize), 'same');
		plot (gain,  char(col(ind)));
		hold on;
	end;
	title ('Gainsol amp from chosen antennas');
	xlabel ('Timeslice'); ylabel ('Arbit units');
	grid;


	gainsrc = [abs(complex(re_gain, im_gain))  srcflux];
	% Show zoomed in timeseries of gain amplitudes and src fluxes.
	figure;
	for ind = 1:5
		subplot (211);
		plot (conv (gainsrc (1:240, 24*ind), 1/winsize*ones (1, winsize), 'same'), ...
			  char (col(ind)));
		hold on;
		subplot (212);
		plot (conv (gainsrc (1:240, 288+ind), 1/winsize*ones (1, winsize), 'same'), ...
			  char (col(ind)));
		hold on;
	end;
	subplot (211);
	xlabel ('Timeslices'); ylabel ('Gain amplitude'); 
	title ('Gain amplitude from chosen antennas (24*ind)');

	subplot (212);
	xlabel ('Timeslices'); ylabel ('Src flux');
	title ('Src flux estimated via LS imaging');

	% Show correlation between model source fluxes and calib. gain solutions
	[r,p] = corrcoef (gainsrc);
	figure; imagesc (r); colorbar;
	title ('Corr coe. between gain amps and LS src flux');
	xlabel ('Ant. number + src number (last 5)');
	ylabel ('Ant. number + src number (last 5)');

	% Store fit parameters 
	coe = zeros (ncal, order+1);
	xcoe = zeros (ncal, 3);
	ycoe = zeros (ncal, 3);
	resid = zeros (ntslices, ncal);
	
	for ind = 1:ncal
		fprintf (1, '--> Operating on source %d\n', ind);
		x = [1:ntslices];    % The x-axis data (time axis)
		% dat = peakfl(ind, :);
		dat0 = conv (peakfl(ind, :), 1/winsize*ones (1, winsize), 'same');

		% Remove outliers
		nout = 0; % number of outliers.
		for jind = 1:5
			dat = dat0 (isnan (dat0) == 0);
			mu = mean (dat);
			sig = std (dat);
			meanmat = repmat (mu, 1, length (dat));
			sigmat = repmat (sig, 1, length (dat));
			out = abs (dat-meanmat) > 3*sigmat;
			fprintf (1, '   Found %d outliers.\n', sum (out));
			% x (out == 1) = [];
			dat (out == 1) = NaN; % [];
			nout = nout + sum (out);
			dat0 (isnan (dat0) == 0) = dat; % Allows keeping outlier pos. safe.
		end;
		fprintf (1, '--> For source %d, total outliers = %d\n', ind, nout);	
		sel = (isnan (dat0) == 0);
		% detrend flux timeseries
		% fit polynomial model of desired order.
		[coe(ind, :) , errest(ind)] = polyfit (x(sel), dat0(sel), order);
		[pfit, del] = polyval (coe(ind, :), x, errest(ind));
		% Estimate errors in polynomial coefficients
		Rinv = inv (errest(ind).R);
		% Covariance matrix of coefficients, see polyfit help.
		coecov = (Rinv*Rinv')*errest(ind).normr^2/errest(ind).df;
		fprintf (1, '--> Polyfit order: %d, coefficients: ',order);
		for jind = 1:order
			fprintf (1, '%f ', coe(ind,jind));
		end;
		fprintf (1, '\n-->Cov matrix: \n');
		disp (coecov);

		% Compute R and R.^2 statistics (TODO).	
		resid (:, ind) = (dat0 - pfit);

		figure (peakplt);
		% Plot confidence bounds
		% Plot the data, the fit, and the confidence bounds
		subplot (3, ncal, ind);
		plot(x, dat0,'-.',...
		     x, pfit,'g-',...
		     x, pfit+2*del,'r*',...
		     x, pfit-2*del,'r*'); 
		title (sprintf ('Src: 3C%s', srclist3CR(callist(ind)).name));
		xlabel('Timeslice');
		ylabel('Extracted flux (Counts)');
		% title('3rd order Polynomial Fit with Confidence Bounds');
		grid on;
		axis tight;
		
		% Plot residues
		subplot (3, ncal, ind + ncal);
		plot (x, resid (:, ind), '-.');
%		L = get (gca, 'XLim'); 
%		set (gca, 'XTick', linspace (L(1), L(2), numxticks));
		title (sprintf ('Resid. fit order: %d', order));
		xlabel ('Timeslice'); ylabel ('Fit residue');
		grid on;
		axis tight;
		
		subplot (3, ncal, ind + 2*ncal);
		redchisq = plotfithist (resid (:, ind), 50, gca);
		title ('resid hist with fit');
		fprintf (1, '\n-->residual hist reduced chisq: %f, residual std: %f\n', redchisq, std (resid(:,ind)));

		% Plot other fit related things.
		figure (posplt);
		subplot (3, ncal, ind);
		plot (x, conv (sigx(ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% plot (x, sigx(ind,:), '.');
		title (sprintf ('Src: 3C%s', srclist3CR(callist(ind)).name));
		xlabel ('Timeslice'); ylabel ('fit major axis (pix)'); 
		axis tight;
%		L = get (gca, 'XLim'); 
%		set (gca, 'XTick', linspace (L(1), L(2), numxticks));

		subplot (3, ncal, ind + ncal);
		plot (x, conv (sigy(ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% plot (x, sigy(ind,:), '.');
		xlabel ('Timeslice'); ylabel ('fit minor axis (pix)'); 
		axis tight;
%		L = get (gca, 'XLim'); 
%		set (gca, 'XTick', linspace (L(1), L(2), numxticks));

		subplot (3, ncal, ind + 2*ncal);
		plot (x, conv (intfl (ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% plot (x, intfl(ind,:), '.');
		xlabel ('Timeslice'); ylabel ('Integrated flux'); 
%		L = get (gca, 'XLim'); 
%		set (gca, 'XTick', linspace (L(1), L(2), numxticks));
		axis tight;
		 
		% Plot positional info.
		% figure (posplt);
		figure;
		[xcoe(ind, :) , xerrest(ind)] = polyfit (x(sel), xpos(ind, sel), 2);
		[pfit, del] = polyval (xcoe(ind, :), x, xerrest(ind));
		% Estimate errors in polynomial coefficients
		Rinv = inv (xerrest(ind).R);
		% Covariance matrix of coefficients, see polyfit help.
		coecov = sqrt (diag(Rinv*Rinv').*xerrest(ind).normr^2./xerrest(ind).df);
		fprintf (1, '--> l-pos Polyfit order: %d, coefficients: ',order);
		for jind = 1:order
			fprintf (1, '%f ', xcoe(ind,jind));
		end;
		fprintf (1, '\n-->Cov matrix: \n');
		disp (coecov);
		% subplot (2, 2*ncal, 2*(ind-1) + 1);

		subplot (221);
		plot(x, xpos(ind,:),'.',...
		     x, pfit,'g-',...
		     x, pfit+2*del,'r:',...
		     x, pfit-2*del,'r:'); 
		title (sprintf ('Src: 3C%s', srclist3CR(callist(ind)).name));
		xlabel('Timeslice');
		ylabel('Fitted l-pos.(Pixels)');
		% title('3rd order Polynomial Fit with Confidence Bounds');
		grid on;
		axis tight;
		
		% subplot (2, 2*ncal, ncal*2+2*(ind-1) + 1);
		subplot (223);
		redchisq = plotfithist (xpos(ind,sel) -pfit(sel), 50, gca);
		fprintf (1, '\n-->residual hist reduced chisq: %f\n', redchisq);
		% hist (xpos(ind,sel) -pfit(sel), 50);
		xlabel ('Timeslice'); ylabel ('Frequency of residue');

		[ycoe(ind, :) , yerrest(ind)] = polyfit (x(sel), ypos(ind, sel), 2);
		[pfit, del] = polyval (ycoe(ind, :), x, yerrest(ind));
		% Estimate errors in polynomial coefficients
		Rinv = inv (yerrest(ind).R);
		% Covariance matrix of coefficients, see polyfit help.
		coecov = (Rinv*Rinv')*yerrest(ind).normr^2/yerrest(ind).df;
		fprintf (1, '--> m-pos Polyfit order: %d, coefficients: ',order);
		for jind = 1:order
			fprintf (1, '%f ', ycoe(ind,jind));
		end;
		fprintf (1, '\n-->Cov matrix: \n');
		disp (coecov);

		% subplot (2, 2*ncal, 2*(ind-1) + 2);
		subplot (222);
		plot(x, ypos (ind, :),'.',...
		     x, pfit,'g-',...
		     x, pfit+2*del,'r:',...
		     x, pfit-2*del,'r:'); 
		title (sprintf ('Src: 3C%s', srclist3CR(callist(ind)).name));
		xlabel('Timeslice');
		ylabel('Fitted m-pos.(Pixels)');
		% title('3rd order Polynomial Fit with Confidence Bounds');
		grid on;
		axis tight;

		% subplot (2, 2*ncal, ncal*2+2*(ind-1) + 2);
		subplot (224);
		redchisq = plotfithist (ypos(ind,sel)-pfit(sel), 50);
		fprintf (1, '\n-->residual hist reduced chisq: %f\n', redchisq);
		% hist (ypos(ind,sel)-pfit(sel), 50);
   end;

   % Plot correlations between residues of various sources.
   % NOTE: Does not work with NaNs in the data. Going for scatter plot instead.
	% Get an intersection of all src. timeseries without a NaN.
	tsnan = zeros (size (resid, 1), 1);
	for ind = 1:ncal
		tsnan = tsnan | isnan (resid (:, ind));
	end;
	% Generate corrcoef.
	resid_nonnan = resid (tsnan == 0, :);
	[r,p] = corrcoef (resid_nonnan);

	figure;
	for ind = 1:ncal
		for jind = ind:ncal
			subplot (ncal, ncal, (ind-1)*ncal+jind);
			scatter (resid (:,ind), resid (:,jind)); 
			axis ([-800 800 -800 800]);
			xlabel ('Resid (count)');
			ylabel ('Resid (count)');
			if ind == jind
				continue;
			else
				subplot (ncal, ncal, (jind-1)*ncal+ind);
				text (0.5, 0.5, sprintf ('r = %f\np = %f', r(ind, jind), p(ind, jind)));
			end;
		end;
	end;
	% Put in correlation coefficients
	for ind = 1:ncal
		subplot (ncal, ncal, ind);
		title (sprintf ('3C%s', srclist3CR(callist(ind)).name));
		subplot (ncal, ncal, 1 + (ind-1)*ncal);
		ylabel (sprintf ('3C%s', srclist3CR(callist(ind)).name));
	end;
	p = mtit ('Scatter plot for monitor source extracted flux timeseries', 'yoff',0.045);
	
	if (fid) fclose (fid); end;
	if (funcal) fclose (funcal); end;
