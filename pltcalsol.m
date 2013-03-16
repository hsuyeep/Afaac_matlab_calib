% Script to plot the calibration solutions available from a .bin file, which
% in turn is generated using the wrcalvis2bin.m script.
% pep/13Jan13
% Added plotting of gain ratios and temporal ratios and differences, and 
% histograms of gain amplitudes and phases. 
% Added plotting of wsf positions Vs. catalog positions.
% pep/20Jan13
%
% Arguments:
%	fname:	file holding gain solutions.
% 	nrecs: Number of timeslices to plot.
% showplt: Bool to control the generation of active plots or not. If 0, only the
% 		   summary timeseries of gain amplitudes and phases are generated.
%
% Returns:
%	re_gain/im_gain: Timeseries of gain solutions.
%	srcflux:	Timeseries of model source fluxes extracted by LSImaging.(NOT WORKING!)

function [srcflux, re_gain, im_gain] = pltcalsol (fname, nrecs, showplt)
	fid = fopen (fname, 'rb');
	Nelem = 288;
	tobs = fread (fid, 1, 'double');	
	fclose (fid);
	fid = fopen (fname, 'rb');
	ind = 1;
	col = {'bo', 'mo', 'ro', 'ko', 'go', 'yo', 'wo', 'co'};

	try
		rec0 = readcalsol (fid);
	catch err
		fprintf (2, 'pltcalsol: Eof reached!\n');
		return;
	end;
	t_first = rec0.tobs;
	if (isempty(rec0.tobs) == 1) 
		disp('End of file reached!'); 
		return;
	end;
	prevgain = 1;
	if (isempty(nrecs) || nrecs < 0)
		% Determine number of records: Crude way, as record size could not be 
		% determined correctly!
		t = whos ('rec0');
		recsize = t.bytes; 
		d = dir (fname);
		nrecs = int32 (d.bytes/t.bytes);
		fprintf (1, '-->Filesize: %d, recsize: %d, nrecs: %d', ... 
				 d.bytes, t.bytes, nrecs);
	end;

	gainsol = complex (rec0.real_gainsol, rec0.imag_gainsol);
	% flagant = rec.flagants; % (abs(gainsol) == 0);
	rem_ants = Nelem - length (rec0.flagants);
	currgain = gainsol; % (flagant == 0);
	gaindiff = currgain - prevgain;
	gainratio = currgain ./ prevgain;

	prevgain = currgain;
	sigman_vec = complex (rec0.real_sigman, rec0.imag_sigman);
	sigman = reshape (sigman_vec, rem_ants, rem_ants);
	caliter = zeros (2, nrecs); % holds number of cal_ext iterations per sol
	stefiter = caliter;

	rec = readcalsol (fid);
	dt = rec.tobs - t_first;


	fseek (fid, 0, 'bof');
	re_gain = zeros (nrecs, Nelem); 
	im_gain = zeros (nrecs, Nelem); 
	srcflux = zeros (nrecs, rec.calsrcs);
	% pause;

	% Figure management
	if (showplt == 1)
		fluxplt = figure;
		gainplt = figure;
		noiseplt = figure;
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

	for ts = 1:nrecs;
		% NOTE: Order of reads is as defined in wrcalsol2bin.m!
		try
			rec = readcalsol (fid);
		catch err
			fprintf (2, 'pltcalsol: Eof reached!\n');
			break; 
		end;
		if (isempty(rec.tobs) == 1) 
			disp('End of file reached!'); break;
		end;
		gainsol = complex (rec.real_gainsol, rec.imag_gainsol);
		rem_ants = Nelem - length (rec.flagants);
		re_gain(ts, :) = rec.real_gainsol;
		im_gain(ts, :) = rec.imag_gainsol;
		caliter (1, ts) = rec.calext_iters;
		caliter (2, ts) = rec.pinv_sol;
		stefiter (1, ts) = rec.stefcal_iters;
		stefiter (2, ts) = rec.stefcal_tol;

		% srcflux (ts, :) = rec.sigmas';
		sigman_vec = complex (rec.real_sigman, rec.imag_sigman);
		sigman = reshape (sigman_vec, rem_ants, rem_ants);

		ind = ind + 1;
		t = whos ('rec');
		fprintf (1, 'Time: %.2f (rec: %03d), recsize: %d, flagants: %03d.', ... 
				 rec.tobs, ts, t.bytes, length(rec.flagants));
		fprintf (1, 'Calextiter:%02d, pinv:%f, stefiter:%02d, tol:%f\n', ... 
  		 rec.calext_iters, rec.pinv_sol, rec.stefcal_iters, rec.stefcal_tol); 

		if (showplt == 1)
			figure (fluxplt);
			for src = 1:rec.calsrcs
				plot (rec.tobs-t_first, rec.sigmas(src), char(col(src)));
				hold on;
			end;
			title (sprintf ('Extracted fluxes of model sources from %s', fname));
			xlabel ('Time (sec from obs. commencement)');
			ylabel ('Flux ratio normalized to CasA flux');
	
			figure (gainplt);
			currgain = gainsol; % (flagant == 0);
			gaindiff = currgain - prevgain;
			gainratio = currgain ./ prevgain;
			prevgain = currgain;
	
			subplot (2,3,1);
			plot (abs (gaindiff), 'ro');
			title (sprintf ('Gain temporal difference abs: Time: %f', rec.tobs));
	
			subplot (2,3,2);
			plot (angle (gaindiff), 'ro');
			title (sprintf ('Gain temporal difference phase: Time: %f', rec.tobs));
	
			subplot (2,3,3);
			plot (abs(currgain), 'bo');
			title ('Dipole gain abs.');
	
			subplot (2,3,4);
			plot (angle(currgain), 'bo');
			title ('Dipole gain phase');
	
			subplot (2,3,5);
			plot (abs(gainratio), 'go');
			title ('Dipole gain ratio abs.');
	
			subplot (2,3,6);
			plot (angle(gainratio), 'go');
			title ('Dipole gain ratio phase');
			% text(.75,1.25, sprintf ('File : %s', fname));
	
			figure (noiseplt);
			subplot (2,2,1);
			plot (rec.tobs-t_first, norm (sigman, 'fro'), 'ro');
			hold on;
			title (sprintf ('Estimated noise norm: Time: %f', rec.tobs));
			% imagesc(abs(sigman)); % estimated noise absolute value
			% colorbar;
			% title (sprintf ('Estimated Noise magnitude. Time: %f', rec.tobs));
	
			subplot (2,2,2);
			plot (rec.tobs-t_first, rec.pinv_sol, 'bo');
			hold on;
			title (sprintf ('Estimated parameter residuals'));
	
			subplot (2,2,3);
			for ind = 1:rec.calsrcs
				plot (ts, rec.thsrc_cat(ind) - rec.thsrc_wsf (ind), char(col(ind)));
				hold on;
			end;
			title (sprintf ('Catalog - WSF zenith angle for %d sources', ...
					rec.calsrcs));
	
			subplot (2,2,4);
			for ind = 1:rec.calsrcs
				plot(ts, rec.phisrc_cat(ind) - rec.phisrc_wsf (ind), char(col(ind)));
				hold on;
			end;
			title (sprintf ('Catalog - WSF azimuth angle for %d sources', ...
					rec.calsrcs));
			% text(.75,1.25, sprintf ('File : %s', fname));
	
			% imagesc (angle(sigman));
			% colorbar;
			% title (sprintf ('Estimated Noise phase. Time: %f', rec.tobs));
			
			% hold on;
			% pause; % (0.1);
		end;
	end;
	fclose (fid);

	% Now generate histograms of each dipoles' gain amplitudes
	gainrehist = hist (re_gain, 100); % 100 bins
	gainimhist = hist (im_gain, 100); % 100 bins
	figure; 
	subplot (2,2,1); imagesc (gainrehist);
	title ('Per antenna gain, real histogram');
	subplot (2,2,2); imagesc (gainimhist);
	title ('Per antenna gain, imag histogram');
	subplot (2,2,3); 
	for ind = 1:size (gainrehist,2);
		plot (gainrehist(:,ind));
		hold on;
	end;
	title ('Per antenna stacked real histogram.');
	subplot (2,2,4); 
	for ind = 1:size (gainimhist,2);
		plot (gainimhist(:,ind));
		hold on;
	end;
	title ('Per antenna stacked imag histogram.');
	% text(.75,1.25, sprintf ('File : %s', fname));

	% Now generate a plot of the gain and phase solutions as a function of time,
	% First for all antennas, and then for a representative antenna from each
	% station.
	amp_ts = hypot (re_gain, im_gain);
	ph_ts = angle (complex (re_gain, im_gain));
	%amp_ts = re_gain;
	%ph_ts = im_gain;

%	figure;
%	subplot (2,1,1);
%	title ('Time series of per-antenna gain solution magnitudes');
%	surf (1:size (amp_ts, 2), 1:size (amp_ts,1), amp_ts);
%
%	subplot (2,1,2);
%	title ('Time series of per-antenna gain solution phases');
%	surf (1:size (amp_ts, 2), 1:size (amp_ts,1), ph_ts);

	figure;
	subplot (2,1,1);
	for ind = 1:6 % 6 stations
		plot (amp_ts(:,24*ind), char(col(ind)));
		title ('Timeseries of gain magnitudes from individual antennas');
		hold on;
	end;
	hold off;

	subplot (2,1,2);
	for ind = 1:6 % 6 stations
		plot (ph_ts(:,24*ind), char(col(ind)));
		title ('Timeseries of gain phases from individual antennas');
		hold on;
	end;
	% text(.75,1.25, sprintf ('File : %s', fname));
	hold off;

	figure;
	subplot (1,2,1)
	hist (caliter(1,:));
	subplot (1,2,2)
	hist (stefiter(1,:));
	% plotyy ([1:nrecs], caliter(1,:), [1:nrecs], stefiter(1,:)) ;
	
