% Script to carry out analysis of source fluxes extracted using func_calfluxrat
% pep/23Jul13
% Arguments:
%	winsize: Size of smoothening tophat kernel to apply.
%   recoff : record offset from which to start processing
%  recs: Number of consecutive timeslices to consider. -1 indicates all.
%   order  : Order of polynomial to fit on flux timeseries for each source.

function func_calfluxratanal (winsize, recoff, recs, order)
	load '../srclist3CR.mat';

	% Night time data
	% Load extracted flux timeseries (generated after running func_calfluxrat.m)
	load ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch_6_convcal_7_el_fftimg.bin_save.mat');
	% Load calibration solution file for solution timeseries.
	 fid = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch_6_convcalsol.bin', 'rb');
	% Load uncalibrated visibility file for pulling out autocorrelations.
	funcal = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch.bin', 'rb');

	% Day time data.
	% Load extracted flux timeseries
%	load ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_7_convcal_6_el_fftimg_fliplr.bin_save.mat');
%%	% Load calibration solution file for solution timeseries.
%	fid = fopen ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_7_convcalsol.bin', 'rb');
%%	% Load uncalibrated visibility file for pulling out autocorrelations.
%	funcal = fopen ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr.bin', 'rb');
	

	% Determine time resolution of data.
	try
		rec0 = readcalsol (fid);
	catch err
		fprintf (2, 'func_calfluxratanal: Eof reached!\n');
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
	[ncal totrecs] = size (fit);
	% Find out skipping used in lightcurve gen.
	fitdt = tobs(2) - tobs(1); 
	% NOTE: Assume calsols are always at the highest resolution!
	recs2skip = fitdt/dt; 
	% if (recs < 0) recs = floor (totrecs*recs2skip - recoff); end;
	if (recs < 0) recs = floor ((tobs(end)-tobs(1))/dt); end;
	fprintf (1, 'Sources: %d, Timeslices (Total/processed): %d/%d\n', ...
			 ncal, totrecs, recs);
	recend = recoff + totrecs;

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

	% Generate required timeseries. Note that NaNs are not plotted, so are good
	% for showing missing data in a timeseries.
	% NOTE: this can grow due to missing data in the light curve!
	peakfl = NaN (ncal, recs); 
	intfl = peakfl;
	xpos = peakfl;
	ypos = xpos;
	sigx = xpos; 
	sigy = xpos;
	pa   = xpos;
	fitresid = xpos;
	% tobs = zeros (recs, 1);
	re_gain = NaN (recs, 288);
	im_gain = NaN (recs, 288);
	acorr   = NaN (recs, 288);
	acorrts = NaN (recs, 1);
	srcflux = NaN (recs, rec.calsrcs);
	sysnoise= NaN (recs, 288);

	ts = 1; % local index for dealing with record offsets.
	ind = recoff;
	while (ind < totrecs-1) % NOTE: ignoring offsets for now.
	% while (ind < recend -1)
	% for ind = recoff:recend - 1

		%%% Read in autocorrelation.
		[acc, acorrts(ts), freq] = readms2float (funcal, -1, -1, 288);

		%%% Read in calibration solution 
		rec = readcalsol (fid);

		% Bring all records to a common time
		alltimes = [acorrts(ts) rec.tobs tobs(ind)];
		if (sum(alltimes - alltimes(1)) ~= 0)
			fprintf (1, '###  Timing mismatch: acorrts:%f, calsol:%f, tobs:%f\n'...
					, acorrts(ts), rec.tobs, tobs(ind));
			
			[tim, mind] = max (alltimes);
			while (tobs(ind) < alltimes(mind))
				ind = ind + 1;
			end;
	
			while (rec.tobs < alltimes(mind))
				rec = readcalsol (fid);
			end;
			while (acorrts(ts) < alltimes (mind))
				ts = ts + 1;
				
				[acc, acorrts(ts), freq] = readms2float (funcal, -1, -1, 288);
			end;
		end;
			
		flagant = zeros (1, 288);
		flagant (rec.flagant) = 1;

		re_gain(ts, :) = real (rec.gainsol);
		im_gain(ts, :) = imag (rec.gainsol);
		srcflux (ts, :) = rec.sigmas';
		sysnoise (ts, flagant == 0) = diag (rec.sigman);
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
		pa   (:,ts) = fitparam (:,6);
		fitresid(:, ts) = fit(:,ind).resid;
		fprintf (1, 'ind: %d flagant: %d, sigman size : %d %d, diag: %d\n',...
				 ind, sum(rec.flagant == 1), size (rec.sigman, 1), ...
				 size(rec.sigman, 2), length (diag(rec.sigman)));

		ts = ts + 1; % increment local time counter
		ind = ind + 1;
	end;

	sel = zeros (length(peakfl), ncal); % Holds outliers per monitor src
	sel0 = ones (length (peakfl), 1);	% Holds intersection of all outliers. 

	% Determine outliers from the timeseries of the first mon. src. 
	for ind = 1:ncal
		fprintf (1, '--> Operating on source %d\n', ind);
		x = [1:length(peakfl)]; % NOTE: Missing values are accounted for here!
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
			dat (out == 1) = NaN; 
			nout = nout + sum (out);
			dat0 (isnan (dat0) == 0) = dat; % Allows keeping outlier pos. safe.
		end;
		fprintf (1, '--> For source %d, total outliers = %d\n', ind, nout);	
		sel (:,ind) = (isnan (dat0) == 0);
		sel0 = sel0 & sel(:,ind); % sel0 = Intersecton of all selections.
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
	[r_acorr, p_acorr] = corrcoef(acorr(sel0==1,:));
	figure; imagesc (r_acorr); colorbar;
	xlabel ('Ant. number'); ylabel ('Antenna number');
	title ('Corr.coe between Autocorrs of all antennas');


	% Show system noise timeseries
	figure;
	for ind = 1:6 % Plot for 6 stations
		gain = conv (abs(sysnoise (:, 24*ind)), 1/winsize*ones (1, winsize),...
					 'same');
		plot (gain,  char(col(ind)));
		hold on;
	end;
	title ('Estimated system noise timeseries from 6 antennas (24*ind)');
	xlabel ('Timeslice number'); ylabel ('System noise (arb. units)');
	grid;
	% Show Correlation between system noise timeseries
	[r_acorr, p_acorr] = corrcoef (abs(sysnoise(sel0==1,:)));
	figure; imagesc (r_acorr); colorbar;
	xlabel ('Ant. number'); ylabel ('Antenna number');
	title ('Corr.coe between system noise of all antennas');


	% Show model source flux timeseries
	figure;
	for ind = 1:rec.calsrcs
		plot (conv (srcflux(:, ind), 1/winsize*ones (1, winsize), 'same'), ...
			  char (col(ind)));
		hold on;
	end;
	title ('Model source LS imaging flux estimate timeseries');
	xlabel ('Timeslices'); ylabel ('Flux (arbit)');
	legend ('3C461 (CasA)', '3C405 (CygA)', '3C144 (TauA)', '3C274 (VirA)', ...
			'Sun');
	grid;


	% Show calib gain sols.
	figure;
	for ind = 1:6 % Plot for 6 stations
		gain = conv (abs(complex (re_gain (:, 24*ind), im_gain (:, 24*ind))),...
					 1/winsize*ones (1, winsize), 'same');
		subplot (211);
		plot (gain,  char(col(ind)));
		hold on;
		subplot (212);
		plot (diff(sign(diff(gain))), char(col(ind)));
		hold on;
	end;
	subplot (211);
	title ('Gainsol amp from chosen antennas');
	xlabel ('Timeslice'); ylabel ('Arbit units');
	grid;
	subplot (212);
	title ('Gainsol amp Derivative from chosen antennas');
	xlabel ('Timeslice'); ylabel ('Arbit units');
	grid;


	% Show single plot containing average gain amplitudes and monitor source
	% light curves.
	gainsrc = [abs(complex(re_gain, im_gain))  srcflux];
	% Average gains over all antennas.
	gainsum = NaN (1, length (peakfl));
	gainsum (sel0==1) = mean ((hypot (re_gain(sel0==1,:),im_gain(sel0==1,:)))');

	figure;
	subplot (length(callist)+1, 1, 1);
	plot (gainsum, 'r*'); grid on; axis tight;title ('Average Gain amplitude');
	for ind = 1:length(callist)
		subplot (length(callist)+1, 1, ind+1);
		% For smoothening
%		dat0 = peakfl (ind, :);
%		dat = peakfl(ind, isnan(peakfl(ind,:)) == 0);
%		dat = conv (dat, 1/winsize*ones (1,winsize), 'same');
%		dat0 (isnan(peakfl(ind,:)) == 0) = dat;
%		plot (dat0); grid on; axis tight;
		plot (peakfl(ind, :), '*-'); grid on; axis tight;
		% title (sprintf ('3C%s', srclist3CR(callist(ind)).name));
		legend (sprintf ('3C%s', srclist3CR(callist(ind)).name));
	end;
	

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
	[r,p] = corrcoef (gainsrc(sel0==1,:));
	figure; imagesc (r); colorbar;
	title ('Corr coe. between gain amps and LS src flux');
	xlabel ('Ant. number + src number (last 5)');
	ylabel ('Ant. number + src number (last 5)');

	% Store fit parameters 
	coe = zeros (ncal, order+1);
	xcoe = zeros (ncal, 3);
	ycoe = zeros (ncal, 3);
	resid = NaN (length (peakfl), ncal);
	
	for ind = 1:ncal
		selsrc = sel (:,ind);
		fprintf (1, '--> Operating on source %d\n', ind);
		x = [1:length(peakfl)]; % NOTE: Missing values are accounted for here!
		dat0 = peakfl(ind,:); dat0 (selsrc == 0) = NaN;
		dat0 = conv (dat0, 1/winsize*ones (1, winsize), 'same');

		% detrend flux timeseries
		% fit polynomial model of desired order.
		[coe(ind, :) , errest(ind)] = polyfit (x(selsrc==1), dat0(selsrc==1), order);
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
		resid (selsrc==1, ind) = (dat0(selsrc==1) - pfit(selsrc==1));

		% Plot confidence bounds
		% Plot the data, the fit, and the confidence bounds
		figure (peakplt);            % Uncomment for table of plots
		subplot (3, ncal, ind);      % Uncomment for table of plots
%		figure;						   % Separate figure per source
%		subplot (131);
		plot(x, dat0,'-*',...
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
		subplot (3, ncal, ind + ncal);% Uncomment for table of plots
%		subplot (132);
		plot (x, resid (:, ind), '*-');
		title (sprintf ('Resid. fit order: %d', order));
		xlabel ('Timeslice'); ylabel ('Fit residue');
		grid on;
		axis tight;
		
		subplot (3, ncal, ind + 2*ncal);% Uncomment for table of plots
		% subplot (133);
		redchisq = plotfithist (resid (:, ind), 50, gca);
		title ('resid hist with fit');
		fprintf (1, '\n-->residual hist reduced chisq: %f, residual std: %f\n', redchisq, std (resid(:,ind)));

		% Plot other fit related things.
		figure (posplt);
		subplot (5, ncal, ind);
		plot (conv (sigx(ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% plot (x, sigx(ind,:), '.');
		title (sprintf ('Src: 3C%s', srclist3CR(callist(ind)).name));
		% xlabel ('Timeslice'); 
		ylabel ('fit major axis (pix)'); 
		axis tight;

		subplot (5, ncal, ind + ncal);
		plot (x, conv (sigy(ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% plot (x, sigy(ind,:), '.');
		% xlabel ('Timeslice'); 
		ylabel ('fit minor axis (pix)'); 
		axis tight;

		subplot (5, ncal, ind + 2*ncal);
		plot (x, conv (pa (ind,:), 1/winsize*ones (1, winsize), 'same'), '.-');
		% xlabel ('Timeslice'); 
		ylabel ('Fit PA'); 
		axis tight;

		subplot (5, ncal, ind + 3*ncal);
		% plot (diff(conv (intfl (ind,:), 1/winsize*ones (1, winsize), 'same')), '.-');
		plot (x, intfl(ind,:), '-*');
		% xlabel ('Timeslice'); 
		ylabel ('Integrated flux'); 
		axis tight;
		 
		subplot (5, ncal, ind + 4*ncal);
		plot (x, fitresid(ind,:), '-*');
		% xlabel ('Timeslice'); 
		ylabel ('2D fit residue'); 
		axis tight;
		% Plot positional info.
		% figure (posplt);
		figure;
		[xcoe(ind, :) , xerrest(ind)] = polyfit (x(selsrc==1), xpos(ind, selsrc==1), 2);
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
		pos = NaN (length (peakfl)); pos (selsrc==1) = xpos(ind, selsrc==1);
		plot(x, pos,'.',...
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
		redchisq = plotfithist (xpos(ind,selsrc==1) -pfit(selsrc==1), 50, gca);
		fprintf (1, '\n-->residual hist reduced chisq: %f\n', redchisq);
		xlabel ('Timeslice'); ylabel ('Frequency of residue');

		[ycoe(ind, :) , yerrest(ind)] = polyfit (x(selsrc==1), ypos(ind, selsrc==1), 2);
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
		pos = NaN (length (peakfl)); pos (selsrc==1) = ypos(ind, selsrc==1);
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
		redchisq = plotfithist (ypos(ind,selsrc==1)-pfit(selsrc==1), 50);
		fprintf (1, '\n-->residual hist reduced chisq: %f\n', redchisq);
		% hist (ypos(ind,sel)-pfit(sel), 50);
   end;

	return;
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
	p = mtit ('Scatter plot for monitor source extracted flux timeseries', 'yoff',0.055);
	
	if (fid) fclose (fid); end;
	if (funcal) fclose (funcal); end;

	function err = residual(par, srcflux, gainsum)
	  % global X Y Zdata;
		model_trend = par(2) + par(1)*gainsum; 
		err = sum((model_trend - srcflux).^2);
