% Script to generate FFT images from calibrated visibilities available in a 
% .bin file.
% pep/18Oct12
% Arguments:
%	fname   :   Name of file containing calibrated or uncalibrated visibilities.
%	ntslices:   Number of timeslices to image. -1 => all timeslices.
%	offset  :   Initial offset in visibility file from which to start imaging.
%               Set to 0 for no offset.
%   skip    :   Number of records to skip, in the interest of faster imaging.
%	posfilename:Name of file containing the local positions of the array 
%               configuration corresponding to the data in fname.
%	mosaic  :   Bool controlling generation of mosaic images. If false,  zenith
%               pointing images (less computation) are generated.
%	caxisrng:   2-element vector with the coloraxis range for display.
%	wr2file :   Bool controlling writing of generated images to file. If true,
%               the input filename with 'img' appended is created and written
%               out using floats. NOTE: If writing to file,images are not shown.
%   elbeam  :   Bool controlling applying an element beam pattern correction.
%               NOTE: Currently only for X-dipole.
%   pol     :   Polarization to choose: 0 => X, 1 => Y, 2 => I(?)
%  radec    :   Flag to control if RA/DEC plane projected images are to be created.
%  accum    :   Flag controlling the generation of accumulated images.
%  Returns:
%	img_l   :   The l-axis of the generated image.
%	img_m   :   The m-axis of the generated image.
%     img   :   The actual image matrix.
%	NOTE :      The values returned are for the last image in the image 
%				timeseries!


function [img_l, img_m, img] =  ... 
    genfftimage (fname, ntslices, offset, skip, posfilename, mosaic, caxisrng,...
				 wr2file, elbeam, pol, radec, accum)
    % genfftimage (fname,ntslices, offset, posfilename, weight, uvcellsize, mosaic, caxisrng, wr2file)
	% Gridding parameters
	gparm.type = 'pillbox';
    gparm.lim  = 0;
    gparm.duv = 0.5;				% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    gparm.Nuv = 500;				% size of gridded visibility matrix
    gparm.uvpad = 512;				% specifies if any padding needs to be added
    gparm.fft  = 1;
	nfacet = 3;
	facetsize = 256;

	if (elbeam == 0)
		elstr = 'noel';
	else
		elstr = 'el';
	end;

	if (wr2file == 0)
		hdl = figure;
		if (radec == 1)
			projimg = figure;
			ra_grid = linspace (0,2*pi, gparm.uvpad)*12/pi; % Convert to hours
			de_grid = linspace (-pi/2,pi/2,gparm.uvpad)*180/pi; % Convert to deg.
			if (accum == 1)
				acc_radecmap = zeros (gparm.uvpad);
				acc_localmap = zeros (gparm.uvpad);
			end;
		end;
	else
		k = strfind (fname, '.bin');
		if (mosaic == 1)
			imgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_mosimg.bin');
		else
			imgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_fftimg.bin');
		end;
		if (exist (imgfname, 'file') == 2)
			fprintf (2, 'Overwriting existing file: %s. Continue? (Ctrl-C to kill)\n', imgfname);
			pause;   % To prevent overwriting already written files!
		end;
		fimg = fopen (imgfname, 'wb');
		if (fimg < 0)
			disp (ferror(fimg));
			return ;
		end;
		disp (['Writing generated images to file :' imgfname]);
	end;

	if (ntslices == -1) 			% Get number of records from the file
		[ntslices, tmin, tmax, dt] = getnrecs (fname); 
	end;

	% Obtain image parameters from visibility file.
	fin = fopen (fname, 'rb');
	[acc, img.tobs, img.freq] = readms2float (fin, offset, -1, 288);
	Nelem = size (acc, 1);
	lambda = 299792458/img.freq; 		% in m.
	% duv = lambda/2;
	
	station = 0;
	% For imaging from a single station
	% acc_tmp = acc;
	% acc = zeros (48, 48);
	% acc = acc_tmp (1:48, 1:48);
	% acc = acc_tmp (49:96, 49:96);
	% acc = acc_tmp (97:144, 97:144);
    dl = (299792458/(img.freq * gparm.uvpad * gparm.duv)); % dimensionless, in dir. cos. units
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * gparm.uvpad / 2;
    % l = [-lmax:dl:lmax-1e-3];
    % m = l;  % Identical resolution and extent along m-axis

    % Local horizon based coordinates of array in ITRF
    load (posfilename, 'posITRF', 'poslocal'); 

	% load 3CR catalog
	load ('srclist3CR.mat');

    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	
	% Generate primary beam for given frequency from full EM station level 
	% simulations.
	if (mosaic == 0)
		% Image current timeslice. Generate a zenith image.
   		[radecmap, img.map, calvis, img.l, img.m] = ... 
		  fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
					gparm, [], [], img.tobs, img.freq, radec);
	else
		% Image current timeslice. Generate a mosaic image.
   		[img.map, img.l, img.m] = ... 
				genmosaic(acc, img.tobs, img.freq, nfacet, facetsize, ...
						uloc(:), vloc(:), posITRF, 0);
	end;
	mask = zeros (size (img.map));
	mask (meshgrid (img.l).^2 + meshgrid(img.m).'.^2 < 0.9) = 1;

	if (elbeam == 1)
		addpath '~/WORK/AARTFAAC/Afaac_matlab_calib/LBA_beam/CS1/';
		if (pol == 1) % Ypol beam.
			elembeam = calculateLBAbeam (img.l, img.m, img.freq, [1:2:96]); 
			fprintf (1, 'Generating X-pol beam.\n');
		else          % Xpol beam. NOTE: Default is X-pol!
			fprintf (1, 'Generating X-pol beam.\n');
			elembeam = calculateLBAbeam (img.l, img.m, img.freq, [2:2:96]); 
		end;
		elembeam = max(max(elembeam(:,:,1))) ./ elembeam (:,:,1);
	else
		elembeam = ones (size (mask));
	end;

	% UNTESTED! Works well only with cellsizes of <10meters.
	% Generate weighting mask on sampled visibilities. Those in a higher density
	% uvcell are weighted down (divided by a larger number).
	weightmask = ones (size (acc));
%	if (weight == 1)
%		uvec = uloc (:); vvec = vloc (:);
%		weightvec = weightmask (:);
%		for uind = 1:max(uvec)/uvcellsize
%			usel = (abs(uvec) > (uind-1)*uvcellsize) & ...
%				   (abs(uvec) < uind*uvcellsize);
%			for vind = 1:max(vvec)/uvcellsize
%				vsel = (abs(vvec) > (vind-1)*uvcellsize) & ...
%					   (abs(vvec)<vind*uvcellsize);
%				sel = usel & vsel;
%				w = sum (sel); % Weight = total number of visibilities in range.
%				weightvec (usel & vsel) = 1/w;
%			end;
%		end;
%		weightmask = reshape (weightvec, [Nelem, Nelem]);
%	end;
	
% figure;
% mesh (weightmask);

	if (skip == 0 || isempty (skip)) skip = 1; end;
	fprintf (1, 'Skipping %d recs.\n', skip);
	ts = 1;
	while (ts < ntslices)
	% for ts = 1:ntslices % Can't accomodate skips within a for loop.

		if (isempty (acc) == false)
			fprintf (1, 'Slice: %3d, Time:%.2f, Freq:%.2f.\n', ... 
					 offset+ts,img.tobs, img.freq); 
			acc_weighted = acc .* weightmask;

			if (mosaic == 0)
				% Image current timeslice. Generate a zenith image.
		   		[radecmap, img.map, calvis, img.l, img.m] = ... 
				  fft_imager_sjw_radec (acc_weighted(:), uloc(:), vloc(:), ... 
							gparm, [], [], img.tobs, img.freq, radec);
				if (accum == 1)
					acc_localmap = acc_localmap + img.map;
					if (radec == 1)
						acc_radecmap = acc_radecmap + radecmap;
					end;
				end;
			else
				% Image current timeslice. Generate a mosaic image.
		   		[img.map, img.l, img.m] = ... 
						genmosaic(acc_weighted, img.tobs, img.freq, nfacet, facetsize, ...
								uloc(:), vloc(:), posITRF, 0);
			end;

			% Comment to not apply element beam.
			% Apply dipole primary beam correction. 
			% NOTE: Multiply with elembeam pulls up the image edges.
			img.map = img.map .* (elembeam.*mask); 
			
			if (wr2file == 0)
				figure (hdl);
				if (accum == 1)
	   	     		imagesc(img.l, img.m, (real(acc_localmap)));
				else
	   	     		imagesc(img.l, img.m, (real(img.map)));
				end;

				if isempty (caxisrng) == 0
					caxis (caxisrng);
				end;
				% [sel, sel_l, sel_m] = overplotcat (img.tobs, srclist3CR, ... 
				% 									 10, hdl, true);
	        	set(gca, 'FontSize', 16);
			
		        title(sprintf('Rec: %d, Time:%s, Freq:%.2f',offset+ts, datestr(mjdsec2datenum(img.tobs)), ...
						img.freq));
		        axis equal
		        axis tight
				% To match orientation with station images
   				set (gca, 'YDir', 'Normal'); 
		   		set (gca, 'XDir', 'Reverse'); 

		        ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
		        xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
		        set(colorbar, 'FontSize', 16);
				if (radec == 1)
					figure (projimg);
					if (accum == 1)
						imagesc (ra_grid, de_grid, acc_radecmap); colorbar;
					else
						imagesc (ra_grid, de_grid, radecmap); colorbar;
					end;
					xlabel ('RA(Hr)'); ylabel ('Dec(deg)');
		        	title(sprintf('Rec: %d, Time:%s, Freq:%.2f',offset+ts, datestr(mjdsec2datenum(img.tobs)), ...
							img.freq));
		        	% axis equal
		        	axis tight
				end;
			else
				% Write image to output file.
				img.pix2laxis = length (img.l);
				img.pix2maxis = length (img.m);
				wrimg2bin (fimg, img);
			end;
	
			% Read in next visibility set.
			for ind = 1:skip
				[acc, img.tobs, img.freq] = readms2float (fin, -1, -1, 288);
			end;
			if (isempty (acc))
				disp ('End of file reached!');
				break;
			end;
			% pause (0.01);
			% pause;
		else
			disp ('File end reached!');
		end;
		ts = ts + skip;
	end;
	% Return out only the last image
	img_l = img.l; img_m = img.m; img = img.map;
	fclose (fin);
	if (wr2file == 1)		
		fclose (fimg);
	end;
