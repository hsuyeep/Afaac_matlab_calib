% Script to extract the flux of listed sources as a function of time.
% pep/09Jan13
function trackflux (fname)
	addpath '~/WORK/AARTFAAC/Afaac_matlab_calib/'
	debug = 3;
	radec = 0;
	duv = 2.5;					% Default, reassigned from freq. of obs. to
								% image just the full Fov (-1<l<1)
	Nuv = 500; %1000            % size of gridded visibility matrix
	uvpad = 512; %1024          % specifies if any padding needs to be added
	nfacet = 3; facetsize = 256;
	load ('../srclist3CR.mat');
	
	% Local horizon based coordinates of array in ITRF
	load ('../poslocal.mat', 'posITRF', 'poslocal'); 
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	
	% location of calibrated visibilities:  Nighttime data
	fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/SB001_LBA_OUTER_4min_1ch_cal.bin';
	[ntimes, tmin, tmax, dt] = getnrecs (fname);
	ntimes = 10; % TEMP, TODO
	fid = fopen (fname, 'rb');
	
	% longitude and latitude of CS002 (AARTFAAC center).
	lon = 6.869837540;                         % longitude of CS002 in degrees
	lat = 52.915122495;                        % latitude of CS002 in degrees 
	
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	tobs_jd = tobs/86400. + 2400000.5; % Convert to JD.
	lambda = 299792458/freq; 		% in m.
	duv = lambda/2;
	dl = (299792458/(freq * uvpad * duv)); % dimensionless, in dir. cos. units
	lmax = dl * uvpad / 2;

	% Figure out how many sources will be extracted
	[calmap, l, m] = genmosaic(acc, tobs, freq, nfacet, facetsize, ...
							 uloc(:), vloc(:), posITRF, 0);
	% Get image statistics, 64-pixel square, 3sigma threshold.
	[dr, sig] = getimagedr (calmap, 32, 3);
	[out,goodfit] = extract (calmap, l, m, tobs, freq, 300, srclist3CR, sig, 5, debug);
	fit = zeros ([size(out), ntimes]);
	 pause;

	for img = 1:ntimes

		% Generate a map from calibrated visibilities
		[calmap, l, m] = genmosaic(acc, tobs, freq, nfacet, facetsize, ...
							 uloc(:), vloc(:), posITRF, 0);
		
		% Get image statistics, 64-pixel square, 3sigma threshold.
		[dr, sig] = getimagedr (calmap, 32, 3);
		
		[fit(:,:, img), goodfit] = ... 
			extract (calmap, l,m, tobs, freq, 300, srclist3CR, sig, 5, debug);
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
		tobs_jd = tobs/86400. + 2400000.5; % Convert to JD.
	end;

	% Common time axis definition for all plots
	t_off = min (squeeze (fit(1,8,:)));
	t_axis = [1:length(squeeze (fit(1,8,:)))];

	% Plot estimated peak fluxes as a function of time for all sources
	for src = 1:goodfit; % size (fit, 1)
		subplot (2, round(size (fit, 1)/2), src)
		plot (t_axis, squeeze (fit (src, 3, :)), '+-'); % Amp
		title (sprintf ('Src: 3C%s', srclist3CR(fit (src, 7,1)).name));
	end;
	
%	% Plot all parameters for each detected source, as a func. of time.
%	for src = 1:1 %size(fit, 1)
%		subplot (2,3,1);
%		plot (t_axis, squeeze (fit (src, 1, :)), '+-'); % xpos
%		ylabel ('xpos (pix)');
%
%		subplot (2,3,2);
%		plot (t_axis, squeeze (fit (src, 2, :)), '+-'); % ypos
%		ylabel ('ypos (pix)');
%
%		subplot (2,3,3);
%		plot (t_axis, squeeze (fit (src, 3, :)), '+-'); % Amp.
%		ylabel ('Amp.');
%
%		subplot (2,3,4);
%		plot (t_axis, squeeze (fit (src, 4, :)), '+-'); % sigx
%		ylabel ('Sigx');
%
%		subplot (2,3,5);
%		plot (t_axis, squeeze (fit (src, 5, :)), '+-'); % sigy
%		ylabel ('Sigy');
%		hold on;
%	end;
%
	save ('trackflux.mat', 'fit', 'goodfit');
