% Script to plot the calibration solutions available from a .bin file, together
% with the visibility corresponding to the solution timeslice.
% Ideally used to examine bad calibration solutions, and the corresponding
% visibilities which generated them.
% pep/27Mar13
% Arguments:
%	fvis :	file name of calibrated or uncalibrated visibilities.
%	fsol : file name of calibration solutions.
% 	nrecs: Number of timeslices to plot.
% 	kbhit: bool indicating whether to pause after every display or not.
%
% Returns:
%	re_gain/im_gain: Timeseries of gain solutions.
%	srcflux:	Timeseries of model source fluxes extracted by LSImaging.(NOT WORKING!)

function pltviscalsol (fvis, fsol, nrecs, kbhit)
	% Open visibilities and calibration solutions.
	fid = fopen (fvis, 'rb');
	fsolid = fopen (fsol, 'rb');
	Nelem = 288;
	ind = 1;
	col = {'bo', 'mo', 'ro', 'ko', 'go', 'yo', 'wo', 'co'};
	waitind = 0;
	radec = 0;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500; %1000                % size of gridded visibility matrix
    uvpad = 512; %1024              % specifies if any padding needs to be added

	% Figure management
	visplt = figure;
	visimg = figure;
	skyimg = figure;
	skycalimg = figure;
	set(0,'Units','pixels') 
	scnsize = get(0,'ScreenSize');

	position = get(visplt,'Position');
	outerpos = get(visplt,'OuterPosition');
	borders = outerpos - position;
	edge = -borders(1)/2;

	pos1 = [edge, scnsize(4) * (1/2), scnsize(3)/2 - edge, scnsize(4)/2];
	set(visimg,'OuterPosition',pos1);

	pos1 = [edge+scnsize(3)/2, scnsize(4) * (1/2), scnsize(3)/2 - edge, ... 
			scnsize(4)/2];
	set(visplt,'OuterPosition',pos1);
		
	pos1 = [edge, 0, scnsize(3)/2 - edge, scnsize(4)/2];
	set(skyimg,'OuterPosition',pos1);

	pos1 = [edge+scnsize(3)/2, 0, scnsize(3)/2 - edge, ... 
			scnsize(4)/2];
	set(skycalimg, 'OuterPosition', pos1);

	% Read in the first calibration solution.
	try
		sol = readcalsol (fsolid);
	catch err
		fprintf (2, 'pltcalsol: Eof reached!\n');
		return;
	end;
	lambda = 299792458/sol.freq; 		% in m.
	duv = lambda/2;
    dl = (299792458/(sol.freq * uvpad * duv)); 
    lmax = dl * uvpad / 2;
    load ('poslocal.mat', 'posITRF', 'poslocal'); 
    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	mask = zeros (uvpad);

	t_first = sol.tobs;
	if (isempty(sol.tobs) == 1) 
		disp('End of file reached!'); 
		return;
	end;
	prevgain = 1;
	if (isempty(nrecs) || nrecs < 0)
		% Determine number of records: Crude way, as record size could not be 
		% determined correctly!
		t = whos ('sol');
		recsize = t.bytes; 
		d = dir (fname);
		nrecs = int32 (d.bytes/t.bytes);
		fprintf (1, '-->Filesize: %d, recsize: %d, nrecs: %d\n', ... 
				 d.bytes, t.bytes, nrecs);
	end;

	% Process each calsol.
	for ts = 2:nrecs
		% Read in a visibility record.
		[acc, tobs_vis, freq] = readms2float (fid, -1, -1, 288);
	
		% align to solution record
		if (tobs_vis ~= sol.tobs) 
			offset = int32 (tobs_vis - sol.tobs);
			fprintf (2, '\n<--Time vis:%.2f  sol: %.2f, offset: %4.1f.', ...
					 tobs_vis, sol.tobs, offset);
			if (offset > 0)
				% Move to the next solution larger than current vis. time.
				for ind = ts:nrecs
					try
						sol = readcalsol (fsolid);
					catch 
						fprintf (2, 'pltviscalsol: EoF on sol file!\n');
						break;
					end;
					% NOTE: We may not have to move all 'offset' records, as
					% we may be missing records ourselves!
					if (sol.tobs >= tobs_vis)
						break;
					end;
				end;
			else
				% Move in visibility file till sol.tobs is reached.
				% fprintf (2, 'Skipping %d recs on file %s.', abs(offset), fvis);
				for ind = 1:abs(offset)
					try
						[acc, tobs_vis, freq] = readms2float (fid, -1, -1, 288);
					catch 
						fprintf (2, 'pltviscalsol: EoF on sol file!\n');
						break;
					end;
					% NOTE: We may not have to move all 'offset' records, as
					% we may be missing records ourselves!
					if (sol.tobs <= tobs_vis)
						 % fprintf (1, '%.2f %.2f\n', tobs_vis, sol.tobs);
						break;
					end;
				end;
			end;
		end;

		% Operate on aligned timeslices.
		% plot suspect ACM
		figure (visimg);
		imagesc (10*log10(abs(acc - diag(diag(acc)))));
		colorbar;
		title (sprintf ('Time: %.2f', tobs_vis));

		% Plot gainsolutions
		figure (visplt);
		subplot (1,2,1);
		plot (abs(sol.gainsol));
		subplot (1,2,2);
		plot (angle(sol.gainsol));
		title (sprintf ('Time: %.2f', sol.tobs));

		% Plot skymaps, uncalibrated and from solution.
		figure (skyimg);
   		[radecmap, img.map, calvis, img.l, img.m] = ... 
		  fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
					duv, Nuv, uvpad, sol.tobs, sol.freq, radec);
		mask (meshgrid (img.l).^2 + meshgrid(img.m).'.^2 < 0.9) = 1;
		imagesc (img.map.*mask);
		colorbar;
		title (sprintf ('Time: %.2f', sol.tobs));
		
%		figure (skycalimg);
%   		[radecmap, img.map, calvis, img.l, img.m] = ... 
%		  fft_imager_sjw_radec (sol.calvis(:), uloc(:), vloc(:), ... 
%					duv, Nuv, uvpad, sol.tobs, sol.freq, radec);
%		mask (meshgrid (img.l).^2 + meshgrid(img.m).'.^2 < 0.9) = 1;
%		imagesc (img.map.*mask);
%		colorbar;
%		title (sprintf ('Time: %.2f', sol.tobs));
		

		if (kbhit == 1)
			pause;
		end;
		% Read in next record. Always do this at the end of the loop!
		try
			sol = readcalsol (fsolid);
		catch err
			fprintf (2, 'pltcalsol: Eof reached!\n');
			return;
		end;
	end;
