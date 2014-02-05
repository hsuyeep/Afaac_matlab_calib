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
	% col = {'bo', 'mo', 'ro', 'ko', 'go', 'yo', 'wo', 'co'};
	% col = {'b.', 'm.', 'r.', 'k.', 'g.', 'y.', 'w.', 'c.'};
	col = {'b*', 'm*', 'r*', 'k*', 'g*', 'y*', 'w*', 'c*'};
	% col = {'b-', 'm-', 'r-', 'k-', 'g-', 'y-', 'w-', 'c-'};

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
		fprintf (1, '-->Filesize: %d, recsize: %d, nrecs: %d\n', ... 
				 d.bytes, t.bytes, nrecs);
	end;

	% flagant = rec.flagant; % (abs(gainsol) == 0);
	rem_ants = Nelem - length (rec0.flagant);
	currgain = rec0.gainsol; % (flagant == 0);
	gaindiff = currgain - prevgain;
	gainratio = currgain ./ prevgain;

	prevgain = currgain;
	caliter = zeros (2, nrecs); % holds number of cal_ext iterations per sol
	stefiter = caliter;
	gainerr = zeros (1, nrecs);
	flagmask = zeros (1, 288);

	rec = readcalsol (fid);
	dt = rec.tobs - t_first;


	fseek (fid, 0, 'bof');
	re_gain = zeros (nrecs, Nelem); 
	im_gain = zeros (nrecs, Nelem); 
	tstamp  = zeros (nrecs);
	srcflux = zeros (nrecs, rec.calsrcs);
	srcpos_th  = zeros (nrecs, rec.calsrcs);
	srcpos_phi  = zeros (nrecs, rec.calsrcs);
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
		rem_ants = Nelem - length (rec.flagant);
		re_gain(ts, :) = real (rec.gainsol);
		im_gain(ts, :) = imag (rec.gainsol);
		tstamp (ts) = rec.tobs - t_first;
		caliter (1, ts) = rec.calext_iters;
		caliter (2, ts) = rec.pinv_sol;
		stefiter (1, ts) = rec.stefcal_iters;
		stefiter (2, ts) = rec.stefcal_tol;

		flag1 = zeros (1, 288); flag1 (rec.flagant) = 1;
		srcsel = (rec.sigmas ~= 0);
		rec0sol = [rec0.gainsol(flag1 == 0); rec0.sigmas(srcsel); rec0.thsrc_wsf(srcsel); rec0.phisrc_wsf(srcsel)];
		recsol = [rec.gainsol(flag1 == 0); rec.sigmas(srcsel); rec.thsrc_wsf(srcsel); rec.phisrc_wsf(srcsel)];
		gainerr (ts) = (100/length(rec0sol))*sum( abs(rec0sol - recsol) ./ abs(rec0sol));

		srcpos_th (ts,:) = rec.thsrc_wsf-rec.thsrc_cat;
		srcpos_phi (ts,:) = rec.phisrc_wsf-rec.phisrc_cat;
		srcflux (ts, :) = rec.sigmas';

		ind = ind + 1;
		t = whos ('rec');
		fprintf (1, 'Time: %.2f (rec: %03d), recsize: %d, flagant: %03d.', ... 
				 rec.tobs, ts, t.bytes, length(rec.flagant));
		fprintf (1, 'Calextiter:%02d, pinv:%f, stefiter:%02d, tol:%f\n', ... 
  		 rec.calext_iters, rec.pinv_sol, rec.stefcal_iters, rec.stefcal_tol); 

		if (showplt == 1)
			figure (fluxplt);
			for src = 1:rec.calsrcs
				plot (rec.tobs-t_first, rec.sigmas(src), char(col(src)));
				hold on;
			end;
			set (gca, 'FontWeight', 'bold');
			title ('Extracted fluxes of model sources.');
			xlabel ('Time (sec from obs. commencement)');
			ylabel ('Flux ratio normalized to CasA flux');
			axis tight; grid on;
			
	
			figure (gainplt);
			currgain = rec.gainsol; % (flagant == 0);
			gaindiff = currgain - prevgain;
			gainratio = currgain ./ prevgain;
	
			subplot (2,3,1);
			% plot (abs (gaindiff), 'ro');
			plot (abs(currgain) - abs(prevgain));
			title (sprintf ('Gain time diff abs: Time: %.2f', rec.tobs));
	
			subplot (2,3,2);
			% plot (angle (gaindiff), 'ro');
			plot (abs(currgain), 'bo');
			title ('Dipole gain abs.');

	
			subplot (2,3,3);
			plot (abs(gainratio), 'go');
			title ('Dipole gain ratio abs.');
	
			subplot (2,3,4);
			plot (180/pi*unwrap(angle (currgain) - angle (prevgain)));
			title (sprintf ('Gain time diff phase(deg): Time: %.2f', rec.tobs));
	
			subplot (2,3,5);
			plot (angle(currgain), 'bo');
			title ('Dipole gain phase(rad)');
	
			subplot (2,3,6);
			plot (angle(gainratio), 'go');
			title ('Dipole gain ratio phase');
			% text(.75,1.25, sprintf ('File : %s', fname));
	
			figure (noiseplt);
			subplot (2,2,1);
			plot (rec.tobs-t_first, norm (rec.sigman, 'fro'), 'ro');
			grid on;
			hold on;
			title (sprintf ('Estimated noise norm: Time: %.2f', rec.tobs));
			% imagesc(abs(sigman)); % estimated noise absolute value
			% colorbar;
			% title (sprintf ('Estimated Noise magnitude. Time: %f', rec.tobs));
	
			subplot (2,2,2);
			plot (rec.tobs-t_first, rec.pinv_sol, 'bo');
			grid on;
			hold on;
			title (sprintf ('Estimated parameter residuals'));
	
			subplot (2,2,3);
			for ind = 1:rec.calsrcs
				plot (rec.tobs-t_first, 180/pi*(rec.thsrc_cat(ind) - rec.thsrc_wsf (ind)), char(col(ind)));
				hold on;
			end;
			grid on;
			title (sprintf ('Cat-WSF ele. resid(deg), %d srcs', ...
					rec.calsrcs));
	
			subplot (2,2,4);
			for ind = 1:rec.calsrcs
				plot(rec.tobs-t_first, 180/pi*(rec.phisrc_cat(ind) - rec.phisrc_wsf (ind)), char(col(ind)));
				hold on;
			end;
			grid on;
			title (sprintf ('Cat-WSF azi. resid(deg), %d srcs', ...
					rec.calsrcs));
			% text(.75,1.25, sprintf ('File : %s', fname));
	
			% imagesc (angle(sigman));
			% colorbar;
			% title (sprintf ('Estimated Noise phase. Time: %f', rec.tobs));
			
			% hold on;
			% pause; % (0.1);
			prevgain = currgain;
		end;
	end;
	fclose (fid);

	% Stamp the plots with filenames
	if (showplt == 1)
		figure (fluxplt);
		tb1 = uicontrol ('style', 'text');
		set (tb1, 'Units', 'characters');
		pos = get (tb1, 'Position');
	    pos(1) = 0; pos (2) = 0; pos(3) = length(fname); pos(4) = 1; 
		set (tb1, 'Position', pos); set (tb1, 'FontSize', 8);
		set (tb1, 'String', fname);

		figure (gainplt);
		tb1 = uicontrol ('style', 'text');
		set (tb1, 'Units', 'characters');
		set (tb1, 'Position', pos); set (tb1, 'FontSize', 8);
		set (tb1, 'String', fname);
		
		figure (noiseplt);
		tb1 = uicontrol ('style', 'text');
		set (tb1, 'Units', 'characters');
		set (tb1, 'Position', pos); set (tb1, 'FontSize', 8);
		set (tb1, 'String', fname);
	end;

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

	% Plot LS imaging extracted source flux 
	figure;
	for ind = 1:length(rec.sigmas)
		% plot (tstamp, srcflux(:,ind), char(col(ind)));
		plot (srcflux(:,ind), char(col(ind)));
		hold on;
	end;
	title ('Estimated A-team flux');
	xlabel ('Timeslices'); ylabel ('Flux (arbit)');
	legend ('3C461 (CasA)', '3C405 (CygA)', '3C144 (TauA)', '3C274 (VirA)', ...
			'Sun');
	grid; 


	figure;
	h(1) = subplot (2,1,1);
	for ind = 1:6 % 6 stations
		plot (tstamp, amp_ts(:,24*ind), char(col(ind)));
		hold on;
	end;
	hold off;
	grid;
	title (sprintf ('Timeseries of complex gains from individual antennas:%s', fname));
	ylabel ('Gain amplitude (arb. unit)');
	% Fix number of ticks on axis
	L = get (gca, 'XLim');
	set (gca, 'XTick', int32(linspace (L(1), L(2), 5)));
	L = get (gca, 'YLim');
	set (gca, 'YTick', int32(linspace (L(1), 25, 5)));
	set (h(1), 'xticklabel', []); % Turn off the x-axis on top plot.

	h(2) = subplot (2,1,2);
	for ind = 1:6 % 6 stations
		plot (tstamp, ph_ts(:,24*ind), char(col(ind)));
		hold on;
	end;
	% text(.75,1.25, sprintf ('File : %s', fname));
	hold off;
	grid;
	ylabel ('Gain phase (rad)');

	% Fix number of ticks on axis
	L = get (gca, 'XLim');
	set (gca, 'XTick', int32(linspace (L(1), L(2), 5)));
	set (h(2), 'YAxisLocation', 'right');
	% L = get (gca, 'YLim');
	% set (gca, 'YTick', (linspace (-pi, pi, 5)));
	% lab = linspace (-pi, pi, 5);
	% set (gca, 'YTickLabel', sprintf ('%5.2f|', lab));
	% set (gca, 'YTick', (linspace (-pi, pi, 5)));
	xlabel (sprintf ('Time offset (in %0.2fsec) from %.2f', dt, t_first));
	% Move amp and phase plots close to each other.
	pos=get(h,'position');
	bottom=pos{2}(2);
	top=pos{1}(2)+pos{1}(4);
	plotspace=top-bottom;
	pos{2}(4)=plotspace/2;
	pos{1}(4)=plotspace/2;
	pos{1}(2)=bottom+plotspace/2;
	set(h(1),'position',pos{1});
	set(h(2),'position',pos{2});
	linkaxes (h, 'x');

	figure;
	subplot (2,2,1)
	plot (tstamp, caliter(1,:), '-ob');
	title ('Major cycle iterations over time');
%	hist (caliter(1,:));
%	title ('Histogram of cal_ext iterations');
	subplot (2,2,2)
	plot (tstamp, caliter(2,:), '-ob');
	title ('Major cycle residuals');
%	hist (stefiter(1,:));
%	title ('Histogram of stefcal iterations');

	subplot (2,2,3)
	plot (tstamp, stefiter(1,:), '-ob');
	title ('Minor cycle iterations over time');
		
	subplot (2,2,4)
	plot (tstamp, stefiter(2,:), '-ob');
	title ('Minor cycle residuals');
	% plotyy ([1:nrecs], caliter(1,:), [1:nrecs], stefiter(1,:)) ;

	figure;
	plot (tstamp, gainerr, '-ob');
	xlabel ('Time slices');
	ylabel ('Gain error (%)');
	title ('Error in gain solutions wrt. first solution');
	L = get (gca, 'XLim');
	set (gca, 'XTick', int32(linspace (L(1), L(2), 5)));
	L = get (gca, 'YLim');
	set (gca, 'YTick', int32(linspace (L(1), 100, 5))); % Not beyond 100%
	grid;

	% Plot residual track of model sources
	modsrc = {'3C461 (CasA)', '3C405 (CygA)', '3C144 (TauA)', '3C274 (VirA)', ...
			'Sun'};
	for ind = 1:rec.calsrcs
		% subplot (2,3,ind);
		% [ax, h1, h2] = plotyy ([1:length(srcpos_th)], 180/pi*srcpos_th(:,ind), [1:length(srcpos_th)], 180/pi*srcpos_phi(:,ind));
		% set (ax(1), 'YLim', [-0.2 0.2]);
		% set (ax(2), 'YLim', [-0.2 0.2]);

		figure;
		subplot (211); plot (180/pi*srcpos_th(:,ind), '-'); grid on; axis tight;
		title (sprintf ('Src %s', modsrc{ind})); ylabel ('Ele. angle (deg)');
		axis ([1, length(srcpos_th), -0.1, 0.1]);

		subplot (212); plot (180/pi*srcpos_phi(:,ind), '-'); grid on;axis tight;
		title (sprintf ('Src %s', modsrc{ind})); ylabel ('Azi. angle (deg)')
	end;
