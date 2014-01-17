% Script to excercise xcallsimg.m
% pep/17Jul13

function fit = func_calfluxrat (fname, offset, ntslices, posfilename, debug)

	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('func_calfluxrat: fid < 0! Quitting.');
		return;
	end;
		
	% Load various meta data
    rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
    rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
    rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	load (posfilename, 'posITRF', 'poslocal');
    rodata.posITRF_fl = posITRF; 	% NOTE: Not removing flagged vis. for now.
    rodata.poslocal = poslocal;
	load srclist3CR.mat

	if (ntslices < 0)
		[ntslices, tmin, tmax, dt] = getnrecs (fname);
	end;

	% NOTE: Current list assumes 3CR catalog!
	% calsrcs = 3C[295, 219, 338, 410.1, 10, 84, 338]
	callist = [200, 133, 237]; %, 218, 4, 54, 237];
	ncal = length (callist);

	fit = zeros (ncal,  ntslices);
	for ind = 1:ntslices
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);

		fprintf (1, '\nRec %4d ', ind);
		[upcals, fittemp] = ...
			xcallsimg (acc, tobs, freq, rodata, ...
						 srclist3CR, callist, true, 0);	
		fit (:,ind) = fittemp.peakfl;	
		% fprintf (1, '%s\n', num2str(fit.peakfl, '%f '));
	end;
	
	% Plot fluxes
	col = {'b*-', 'm*-', 'r*-', 'k*-', 'g*-', 'y*-', 'w*-', 'c*-'};
	for ind = 1:ncal
		leg_str {ind} = srclist3CR(callist(ind)).name;
	end;

	figure;
	subplot (121);
	for ind = 1:ncal
		plot (fit (ind, :), char(col(ind)));
		hold on;
	end;
	title ('Extract peak flux of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);
	subplot (122);
	for ind = 1:ncal
		plot (fit(ind, :)./fit(1,:), char(col(ind)));
		% plot (resid(ind,:), char(col(ind)));
		hold on;
	end;
	title ('Extract peak flux ratio of calibrators'); 
	% title ('resid flux ratio of calibrators'); 
	xlabel ('Timeslices'); ylabel ('Counts');
	legend (leg_str);

	% Plot histograms. 
	figure;
	for ind = 1:ncal
		subplot (2,3,ind);
		hist (fit(ind,:));
		title (sprintf ('Peak flux (src:%s)', ...
				srclist3CR(callist(ind)).name));

		subplot (2,3,ind+ncal);
		hist (fit(ind,:)./fit(1,:));
		title (sprintf ('Peak flux ratio to %s (src:%s)', ...
				srclist3CR(callist(1)).name, srclist3CR(callist(ind)).name));
	end;
