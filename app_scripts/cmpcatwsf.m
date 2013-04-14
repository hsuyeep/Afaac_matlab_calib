% Script to generate plots of catalog Vs. WSF positions of model sources, 
% based on the solutions generated and stored by wrcalvis2bin.m.
% pep/080Apr13
%
% Arguments:
%  fname : Name of calibration solution file.
%  nrecs : Number of timeslices to plot.
%
% Returns:
%  ele_diff, az_diff: Ele., azi differences between catalog and wsf positions.

function cmpcatwsf (fname, nrecs)
	fid = fopen (fname, 'rb');
	srcname = {'3C461', '3C405', '3C144', '3C274', 'Sun'};
	% To rotate coordinates in ITRF to the plane of CS002
	rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ... 
   			   0.9928230000, -0.0954190000, 0.0720990000; ...
   	           0.0000330000,  0.6030780000, 0.7976820000];

	try
		rec0 = readcalsol (fid);
	catch err
		fprintf (2, 'pltcalsol: Eof reached!\n');
		return;
	end;

	
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
	
	% Generate data structures
	t_samp = zeros (1, nrecs);

	% Using NaNs to root out sources which are below the horizon. 
	% Did not work as comparison of 2 NaNs returns a false! 
	srccat_azi= zeros (rec0.calsrcs, nrecs); 
	srccat_el = srccat_azi;
	srcwsf_azi= srccat_azi;
	srcwsf_el = srccat_azi;
	% src3 = zeros (10, nrecs);

	rec = readcalsol (fid);
	dt = int32(rec.tobs - rec0.tobs);
	fseek (fid, 0, 'bof');

	% Main loop
	for ts=1:nrecs
		try
			rec = readcalsol (fid);
		catch err
			fprintf (2, 'cmpcatwsf: Eof reached!\n');
			return;
		end;
		t_samp (ts) = rec.tobs - rec0.tobs;
		
		% Convert WSF and catalog th/phi (Equatorial coordinates in radians), 
		% to ITRF.
		% 1. First convert to cartesian coords. from spherical.
		srcposwsf_x = cos (rec.thsrc_wsf) .* cos (rec.phisrc_wsf);
		srcposwsf_y = sin (rec.phisrc_wsf) .* cos (rec.thsrc_wsf);
		srcposwsf_z = sin (rec.thsrc_wsf);
		
		srcposcat_x = cos (rec.thsrc_cat) .* cos (rec.phisrc_cat);
		srcposcat_y = sin (rec.phisrc_cat) .* cos (rec.thsrc_wsf);
		srcposcat_z = sin (rec.thsrc_cat);
		
		% Collect xyz coordinates for individual sources together
		srcwsf_xyz = [srcposwsf_x srcposwsf_y srcposwsf_z];
		srccat_xyz = [srcposcat_x srcposcat_y srcposcat_z];

		% 2. Rotate to CS002
		srcwsf_rot = srcwsf_xyz * rotmat;
		srccat_rot = srccat_xyz * rotmat;
		
		% 3. Calculate ele. and azi.
		sel = rec.sigmas > 0.01*rec.sigmas(1); % Take sources > 1% of CasA flux.
		srccat_el (sel,ts) = asin (srccat_rot (sel,3));
		srccat_azi (sel,ts) = atan2 (srccat_rot (sel,1), srccat_rot (sel,2));
		
		srcwsf_el (sel,ts) = asin (srcwsf_rot (sel,3));
		srcwsf_azi (sel,ts) = atan2 (srcwsf_rot (sel,1), srcwsf_rot (sel,2));
%		for ind=0:4
%			src3 (2*ind+1, ts) = rec.thsrc_wsf(ind+1);
%			src3 (2*ind+2, ts) = rec.phisrc_wsf(ind+1);
%		end;

	end;

	figure;
	% axis ([-3*180/pi 1.4*180/pi 0.3*180/pi 1*180/pi]);
	for src=1:size (srccat_azi, 1)-1 % Don't plot the Sun
		sel = srccat_azi (src,:) ~= 0; % NOTE: This comparison fails for NaNs!
%		if (srccat_azi (src, sel) < -pi)
%			srccat_azi (src, sel) = srccat_azi (src, sel) + pi;
%		end;
%		if (srccat_el (src, sel) < -pi)
%			srccat_el (src, sel) = srccat_el (src, sel) + pi;
%		end;
%
%		if (srcwsf_azi (src, sel) < -pi)
%			srcwsf_azi (src, sel) = srcwsf_azi (src, sel) + pi;
%		end;
%		if (srcwsf_el (src, sel) < -pi)
%			srcwsf_el (src, sel) = srcwsf_el (src, sel) + pi;
%		end;
		plot (srccat_azi(src, sel)*180/pi, ... 
			  srccat_el(src, sel)*180/pi, '.b'); 
		%	  'LineWidth', 3);
		hold on
		plot (srcwsf_azi(src, sel)*180/pi, ...
			  srcwsf_el(src, sel)*180/pi,  '.r');
		%	  'LineWidth', 3);
	end;
	set(gca, 'FontSize', 16);
	title (sprintf ('Predicted Vs. observed track. %d recs at %f dt.', ts, dt));
	xlabel ('Azimuth (deg)');
	ylabel ('Elevation (deg)');
	tb1 = uicontrol ('style', 'text');
	set (tb1, 'Units', 'characters');
	pos = get (tb1, 'Position');
	pos(1) = 0; pos (2) = 0; pos(3) = length(fname); pos(4) = 1; 
	set (tb1, 'Position', pos); set (tb1, 'FontSize', 8);

	% PLot histograms of differences
	for src=1:size (srccat_azi, 1)-1 % Don't plot the Sun
		sel = srccat_azi (src,:) ~= 0; % NOTE: This comparison fails for NaNs!
		if (sum(sel) == 0)
			continue;
		end;
		figure;
		% subplot (1,2,1);
		dist = sqrt ( ...
				((srcwsf_azi (src, sel) - srccat_azi (src, sel))*180/pi).^2 ...
			 +	((srcwsf_el (src, sel) - srccat_el (src, sel))*180/pi).^2);
		[m, v, sel1] = robustmean (dist, 5);

		subplot (2,1,1);
		plot (t_samp(sel), dist, '.');
		title (sprintf ('Deviation between catalog and WSF positions. src: %s, mean:%.2f, sig:%.2f', char(srcname(src)), m, v));

		subplot (2,1,2);
		hist (dist, 100);
	end;
