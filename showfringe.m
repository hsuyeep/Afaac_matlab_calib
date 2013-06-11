% Script to plot the fringes from selected baselines. Specify the filename, 
% and a baseline range in lambda units.
% pep/24May13

function showfringe (fname, posfname, blinerange, fftsize) 

	if (isempty (fname))
		error ('showfringe: specify filename');
	end;

	[nrecs, tmin, tmax, dt] = getnrecs (fname);
	colp = {'b.', 'm.', 'r.', 'k.', 'g.', 'y.', 'w.', 'c.'};
	coll = {'b-', 'm-', 'r-', 'k-', 'g-', 'y-', 'w-', 'c-'};

	% load position information
	load (posfname, 'poslocal');
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
	
	fid = fopen (fname, 'rb');
	[acc, tfirst, freq] = readms2float (fid, 1, -1, 288);
	lambda = 299798456/freq; % In m.
	uloc = uloc ./ lambda; vloc = vloc ./ lambda; % In lambda units now.
	uvdist = sqrt (uloc.^2+vloc.^2);
	% sel = zeros (size (acc)); 
	sel = (triu(uvdist) > blinerange(1)); % && uvdist < blinerange (2));
	% sel (1,288) = 1;

	fprintf (2, 'showfringe: Available lambda range:[%.2f %.2f]\n', ...
			 min(min(uvdist)), max(max(uvdist)));
	if (sum(sum(sel)) == 0) % No baselines selected
		fprintf (2, 'No baselines selected by criteria! Selecting longest 2 percent of the baselines.\n');
		sel = (triu(uvdist) > 0.98*max(max(uvdist)));
	end;
	fprintf (2, 'showfringe: Number of selected baselines: %d\n', ...
			 sum(sum(sel)));

	
	selblines = zeros (sum(sum(sel)), fftsize);
	window = zeros (1, fftsize);
	for ind = fftsize/2:fftsize-1
		window(ind) = (fftsize-ind+1)/(fftsize/2);
		window(fftsize-ind) = window(ind); % Reflect the window.
	end;
	
	tobsblk = zeros (1, fftsize);
	% Generate fringe plots for selected baselines
	fringeplt = figure();
	spectplt = figure ();
	for ts = 1:fftsize:nrecs
		for samp = 1:fftsize
			[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
			if (isempty (tobs)) 
				break; 
			end;

			tobsblk (samp) = tobs;
			selblines(:, samp)  = acc(sel);
		end;

		% Plot the fringe as a function of time
		figure (fringeplt);
		for ind = 1:sum(sum(sel));
			subplot (2,2,1);
			plot (tobsblk-tfirst, abs (selblines(ind,:)), char(coll(ind)));
			hold on
			subplot (2,2,2);
			plot (tobsblk-tfirst, angle (selblines(ind,:)), char(colp(ind)));
			hold on;
			subplot (2,2,3);
			plot (tobsblk-tfirst, real (selblines(ind,:)), char(coll(ind)));
			hold on;
			subplot (2,2,4);
			plot (tobsblk-tfirst, imag (selblines(ind,:)), char(coll(ind)));
			hold on;
		end;

		% Plot the fringe as a function of temporal frequency
		figure (spectplt);
		for ind = 1:sum(sum(sel));
			% subplot (2,1,1);
			ps = fft (selblines (ind, :).*window);
			plot (abs (ps(2:fftsize/2))); % , char(col(ind)));
			hold on;
		end;
		
		pause;
	end;
