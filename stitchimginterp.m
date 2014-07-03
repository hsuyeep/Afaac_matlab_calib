% Script to generate a stitched image from snapshots covering a range of RA 
% positions.
% pep/30May14
% NOTE: Currently works only on a single frequency bin.
% Arguments: 
%  fname: image bin file (as generated by genfftimage
%  skip : Number of images (hence time slices) to skip while creating the stitched 
%		  image.
%		   -1 => Program chooses a time skip factor resulting in a 1 pixel movement
%				 between image fragments.
%		   -2 => Program chooses a time skip factor such that individual fragments 
%			     (corresponding to the field of view fov) are adjacent to each other.
% offset: Timeslice offset from which to start imaging.
% nrecs : The number of timeslices to image. -1 implies all records.
%  fov  : Fraction of field of view to show (0-1), symmetrically about zenith.
%  deb  : Bool debug flag for showing intermediate outputs.

% Returns:
%  st_tim : Timestamps of centers of stitched fragments.
%  st_ra  : RA of timestamps
%  st_dec : Dec of timestatmps
%  st_img : The stitched image.

function [st_tim, st_ra, st_dec, st_img, gridimg] = stitchimginterp (fname, skip, offset, nrecs, fov, deb)

	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('stitchimages: fid < 0! Quitting.');
		return;
	end;

	% Move to correct location within file
	% TODO

	% Read an image to determine relevant imaging parameters.
	img0 = readimg2bin (fid, 0);
	img  = readimg2bin (fid, 0);
	dt = img.tobs - img0.tobs;

	% Determine the number of records provided by the user.
	if (nrecs < 0)
		% Find out the number of records to use in stitching
		recsize = 8+ ... % tobs
				  8+ ... % freq
				  4+ ... % pix2laxis
				  4*img0.pix2laxis + ... % laxis
				  4+...  % pix2maxis
				  4*img0.pix2maxis + ... % maxis
			      4*img0.pix2laxis*img0.pix2maxis; % img contents.
		st = fseek (fid, -recsize, 'eof');
		img_end = readimg2bin (fid, 0);
		nrecs = floor ((img_end.tobs - img0.tobs)/dt);	
		fprintf ('--> Total number of records in file: %d\n', nrecs);
		fclose (fid); 
		fid = fopen (fname, 'rb');
	end;

	% Determine the resolution of a pixel, in deg.
	% NOTE: Only works for small angle approx., so for close-to-zenith pixels.
	% One pixel = pix2tim secs. of time. 12 hrs=1 FoV.
	pix2tim = ((12*60*60)/img0.pix2laxis); 
	recs2pix = pix2tim/dt; % One record = recs2pix pixels.

	% Generate RA/DEC grid. Ensure enough spatial resolution for zenith regions
	% rarng is in rad.
	radecres = pix2tim*(2*pi/86400); % Grid resolution in radians.
	rarng  = [0:radecres:2*pi];
	decrng = [0:radecres:pi];
	% [ra, dec] = meshgrid (rarng, decrng);
	gridimg = zeros (length(decrng), length(rarng));
	
	% Figure out number of pixels around zenith to extract from each snapshot image
	% for creating the stitched image.
	lmin = 1; lmax = img.pix2laxis;  % Default is the full FoV.
	mmin = 1; mmax = img.pix2laxis; 
	% Get dl and dm from the central part of the image.
	dl = img.l(floor(img.pix2laxis/2)) - img.l(floor(img.pix2laxis/2)-1);
	dm = img.m(floor(img.pix2maxis/2)) - img.m(floor(img.pix2maxis/2)-1);
	if (fov < 1)
		lpix2fov = floor ((fov*img.pix2laxis)); 
		mpix2fov = floor ((fov*img.pix2maxis)); 
		lmin = floor(img.pix2laxis/2-lpix2fov/2); % In pixel index units.
		lmax = lmin + lpix2fov;
		mmin = floor(img.pix2maxis/2-mpix2fov/2);
		mmax = mmin + mpix2fov;
	end;

	if (skip == 0) skip = 1; end;
	 
	if (skip == -1) % User requests to skip recs leading to a 1pix movement of frag
		skip = floor(recs2pix);
		fprintf (1, '--> Skipping %d time records to move image by %.3f pixel.\n',...
				 skip, skip/recs2pix);
	elseif (skip == -2) % User requests to skip recs such that individual fovs abut.
		skip = floor ((lmax-lmin+1)*recs2pix); % Integer number of records to skip,
											   % in order to move by (lmax-lmin) pix.
	else
		fprintf (1, '--> User skip value: %d recs/%.3f pixels\n', ...
				 skip, skip/recs2pix);
		skip = floor (floor(skip/recs2pix)*recs2pix);
		fprintf (1, '--> Rounding to integer pixels and recs: %d pix/%d recs\n', ...
				 skip/recs2pix, skip);
	end;

	% Abort if insufficient number of timeslices available.
	if (skip > nrecs)
		error ('stitchimages: Insufficient recs for generating fragments!');
	end;

	if (deb > 0)
		debhdl1 = figure;
		debhdl2 = figure;
	end;


	% Determine the size of the output stitched image.
	% NOTE: This assumes an l/m pixel is the same as an RA/DEC pixel (small FoV)
	nfrags = floor (nrecs/skip)+1;
	nfragcols = lmax-lmin+1; % Number of columns in a single fragment.
	nstrows = mmax-mmin+1; % Number of rows in the stitched image (Pixel units). 
						   % m is the north-south (so declination) axis, and 
						   % hence remains the same for all snapshots.
	nstcols = nfragcols + ceil(skip/recs2pix)*nfrags; % Number of columns in the stitched image (Pixels). 

	% Determine the RA/DEC of the center of this snapshot.
	[rabeg, decbeg] = lmtoradec (0, 0, JulianDay (mjdsec2datenum (img.tobs)), 6.869837540, 52.915122495);

	raend = rabeg; % Placeholder

	% Create datastructures
	st_img = zeros (nstrows, nstcols);
	st_ra = zeros (nstrows, nstcols);
	st_dec = zeros (nstrows, nstcols);
	% Temporary placeholder for ra/dec axis of a fragment.
	frag_ra  = zeros (nstrows, nfragcols);
	frag_dec = zeros (nstrows, nfragcols);
	frag_tim = zeros (1, nstcols);
	frag_l = img.l([lmin:lmax]);
	frag_m = img.m([mmin:mmax]);

	fprintf (1, '\nInput snapshot details:\n');
	fprintf (1, '\tSize: %dx%d pix\n\tFreq:%f Hz\n\ttobs:%s\n', img.pix2laxis, img.pix2maxis, img.freq, datestr (mjdsec2datenum(img.tobs)));
	fprintf (1, '\tpixel: %.3fsecs, %fdeg, %f recs\n', pix2tim, pix2tim*(180/(12*60*60)), recs2pix);
	fprintf (1, 'Fragment details:\n');
	fprintf (1, '\tFragment size: %dx%d pix\n', nstrows, nfragcols);
	fprintf (1, '\tFragment skip: %d recs, %.4f pix\n', skip, skip/recs2pix);
	fprintf (1, '\tNumber of frags in file: %d\n', nfrags);
	fprintf (1, 'Stitched image:\n');
	fprintf (1, '\tSize: %dx%d pix.\n\tRA range: %f-%f\n', nstrows, nstcols, rabeg, raend);


	% For every image in timeseries
	fragind = 1;
	for ind = 1:skip:nrecs
		% Extract out fragment from image
		frag = img.map (lmin:lmax, mmin:mmax);
		% Figure out RA/Dec of the fragment center
		frag_tim (fragind) = img.tobs;
		[frag_ra, frag_dec] = ...
  lmtoradec (frag_l, frag_m, JulianDay (mjdsec2datenum (img.tobs)), 6.869837540, 52.915122495);

		% Determine the starting pixel in the stitched image at which to place this 
		% fragment.
		fragpix = floor((img.tobs - img0.tobs)/pix2tim) + 1;
		fragcolrng = (fragpix:fragpix+nfragcols-1);


		% Find the range of RA/DEC of this patch, in order to interpolate onto 
		% the regular RA/DEC grid.
		frag_ra_rng = floor (min(min(frag_ra))/radecres)+1:floor(max(max(frag_ra))/radecres)+1;
		frag_dec_rng= floor (min(min(frag_dec))/radecres):floor (max(max(frag_dec))/radecres);

   		radecimage = TriScatteredInterp (double(frag_dec(:)),double(frag_ra(:)), ...
										double(real(frag(:))));
		[ra_interp, dec_interp] = meshgrid (double(rarng(frag_ra_rng)), double (decrng(frag_dec_rng)));
		interp_frag = radecimage (double (dec_interp(:)), double (ra_interp(:)));
		interp_frag (isnan(interp_frag)) = 0;
		interp_frag = reshape (interp_frag, [length(frag_dec_rng) length(frag_ra_rng)]);
		gridimg (frag_dec_rng, frag_ra_rng) = gridimg (frag_dec_rng, frag_ra_rng) + interp_frag;
		% Place in stitched image. We assume RA and Dec increments are the same 
		% as an l/m increment.
		st_img (:, fragcolrng) = st_img(:, fragcolrng) + frag;
		st_ra (:, fragcolrng) = frag_ra; 
		st_dec(:, fragcolrng) = frag_dec;
		st_tim = frag_tim;
		if (deb > 0)
			fprintf (1, 'Fragment tim: %f\n', img.tobs);
			fprintf (1, 'Fragment col: (%d:%d).\n', fragpix,  fragpix+nfragcols-1);
%      		[alpha, delta] = lmtoradec (frag_l, frag_m, JulianDay (mjdsec2datenum (img.tobs)), 6.869837540, 52.915122495);
%      		radecimage = TriScatteredInterp (double(alpha (:)), double(delta (:)), ... 
%										double(real(frag(:))));
%			[ragrid, decgrid] = meshgrid (double(alpha(:)), double(delta(:)));
%			radecmap = radecimage (ragrid, decgrid);
			% figure (debhdl1); imagesc (st_img);
			figure (debhdl2); 
			imagesc (rarng*12/pi, decrng*180/pi-90, gridimg); colorbar;
			xlabel ('RA (Hrs)'); ylabel ('Dec (deg)');
			drawnow ();
		end;
		img = readimg2bin (fid, skip);
		fragind = fragind + 1;
		fprintf (2, '.');
		% if (mod (fragind, 10) == 9) fprintf (1, '.'); end;
	end;