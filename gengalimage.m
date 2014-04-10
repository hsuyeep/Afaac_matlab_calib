% Script to incorporate galactic emission from short baselines into calibrated
% visibilities and created images including diffuse emission.
% Created mostly to generate Ralph's PR movie.
% pep/30Aug13

% Arguments:
%  fcal  : File name of calibrated visibilities
%  fsol  : The calibration solution file corresponding to the calib. vis.
% ntslices: Number of timeslices to image
% offset : Time offset from where to generate images
%  skip  : Number of records to skip, for periodic imaging.
% posfilename: The antenna position matrix.
% mosaic : Turn on mosaicing.
% caxisrng: color range for live display.
% wr2file: Write generated images to file or not.
% elbeam : Bool controlling the application of an element beam pattern correction.
% pol    : Beam polarization to choose from: 0 => X, 1 =>Y.

function gengalimages (fnamecal, fnamesol, ntslices, offset, skip, ...
					   posfilename, mosaic, caxisrng, wr2file, elbeam, pol)

	% Open both calib. visibilites and solution files.
	fcal = fopen (fnamecal, 'rb');
	if (fcal < 0)
		fprintf (2, 'Unable to open file %s\n', fnamecal); return;
	end;
	fsol = fopen (fnamesol, 'rb');
	if (fsol < 0)
		fprintf (2, 'Unable to open file %s\n', fnamesol); return;
	end;

	% Compare timestamps
	[acc, img.tobs, img.freq] = readms2float (fcal, offset, -1, 288);
	for ind = 1:offset
		rec = readcalsol (fsol);
	end;
	if (img.tobs ~= rec.tobs)
		fprintf ('Timestamps dont match!\n');
	end;

	
	% Imaging parameters
	gparm.type = 'pillbox';
	gparm.duv = 2.5; 			% Default, reassigned from freq. of obs. to
								% image just the full Fov (-1<l<1)
	gparm.Nuv = 1000;
	gparm.uvpad = 1024; 
	gparm.fft = 1;
	nfacet = 3;
	facetsize = 256;

	Nelem = size (acc, 1);
	lambda = 299792458/img.freq; 	% in m.
    dl = (299792458/(img.freq * uvpad*duv)); % dimensionless, in dir. cos. units
    lmax = dl * uvpad / 2;
    load (posfilename, 'posITRF', 'poslocal'); 
	load ('srclist3CR.mat');
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

	if (ntslices == -1) 			% Get number of records from the file
		[ntslices, tmin, tmax, dt] = getnrecs (fnamecal); 
	end;

	% Generate antenna mask
	antmask = zeros (size (acc));
    posmask = zeros (size (posITRF));
    rem_ants = length(acc) - length(rec.flagant);
    for ind = 1:length(rec.flagant)
        antmask (rec.flagant(ind), :) = 1; antmask (:,rec.flagant(ind)) = 1;
        posmask (rec.flagant(ind), :) = 1;
    end;
	posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]);
	[uloc_fl, vloc_fl] = gen_flagged_uvloc (uloc, vloc, rec.flagant);

	% Generate filename to write into.
	if (elbeam == 0)
		elstr = 'noel';
	else
		elstr = 'el';
	end;

	if (wr2file == 0)
		hdl = figure;
	else
		k = strfind (fnamecal, '.bin');
		if (mosaic == 1)
			imgfname = sprintf ('%s_%d_%s%s', fnamecal(1:k-1), offset, elstr, '_mosimg_gal.bin');
		else
			imgfname = sprintf ('%s_%d_%s%s', fnamecal(1:k-1), offset, elstr, '_fftimg_gal.bin');
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

	% Carry out a single image cycle to get things working...
	% Get rid of flagged antennas
   	acc1 = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
	if (mosaic == 0)
		% Image current timeslice. Generate a zenith image.
	 	[radecmap, img.map, calvis, img.l, img.m] = ... 
			  fft_imager_sjw_radec (acc1(:), uloc_fl(:), ...
					vloc_fl(:), gparm, [], [], img.tobs, img.freq, 0);
	else
		% Image current timeslice. Generate a mosaic image.
   		[img.map, img.l, img.m] = ... 
			genmosaic (acc1, img.tobs, img.freq, nfacet, ...
						facetsize, uloc_fl(:), vloc_fl(:), posITRF_fl, 0);
	end;

%	[radecmap, img.map, calvis, img.l, img.m] = ... 
%				  fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
%							duv, Nuv, uvpad, img.tobs, img.freq, 0);

	% Generate image mask
    imgmask = NaN(length (img.l));
    imgmask (meshgrid(img.l).^2 + meshgrid(img.m).'.^2 < 1) = 1;

	% Generate inner taper
	parm.type = 'Gaussian';
    parm.minlambda = 1;    % NOTE: Units of lambda. 
    parm.maxmeters = 350;  % NOTE: Units of meters.
    parm.pa(1) = 3;       % NOTE: Units of lambda. 
    parm.pa(2) = parm.pa(1);% Inner taper sigx/sigy
    parm.pa(3) = -1;       % NOTE: Units of lambda.
    parm.pa(4) = parm.pa(3);% Outer taper sigx/sigy
    [intapfull, outtap, den, mask, uvdist] =  ...
                taper (posITRF, parm, -1, img.freq, 0);
	intap = reshape (intapfull(antmask ~= 1), [rem_ants, rem_ants]);

	% Generate the element beam pattern if requested for
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
		elembeam = ones (size (imgmask));
	end;

	if (skip == 0 || isempty (skip)) skip = 1; end;
	fprintf (1, 'Skipping %d recs.\n', skip);

	ts = 1;
	while (offset+ts < ntslices)
		if (isempty (acc) == false)
			fprintf (1, 'Slice: %3d, Time:%.2f, Freq:%.2f.\n', ... 
					 offset+ts,img.tobs, img.freq); 
			
			currremants = length (acc)-length (rec.flagant);
			if (currremants ~= rem_ants)
				% Generate new antenna mask
				antmask = zeros (size (acc));
			    posmask = zeros (size (posITRF));
			    rem_ants = length(acc) - length(rec.flagant);
			    for ind = 1:length(rec.flagant)
			        antmask (rec.flagant(ind), :)=1; 
					antmask (:, rec.flagant(ind)) = 1;
			        posmask (rec.flagant(ind), :) = 1;
			    end;
				posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]);
				[uloc_fl, vloc_fl] = gen_flagged_uvloc (uloc, vloc, rec.flagant);
				intap = reshape (intapfull(antmask ~= 1), [rem_ants, rem_ants]);
			end;

			% Get rid of flagged antennas
    		acc = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);

			% Add in the short baselines.
			acc_weighted = acc +  rec.sigman .* intap;

			if (mosaic == 0)
				% Image current timeslice. Generate a zenith image.
		   		[radecmap, img.map, calvis, img.l, img.m] = ... 
					  fft_imager_sjw_radec (acc_weighted(:), uloc_fl(:), ...
							vloc_fl(:), gparm, [], [], img.tobs, img.freq, 0);
			else
				% Image current timeslice. Generate a mosaic image.
   				[img.map, img.l, img.m] = ... 
					genmosaic(acc_weighted, img.tobs, img.freq, nfacet, ...
								facetsize, uloc_fl(:), vloc_fl(:), posITRF_fl, 0);
			end;

			% Comment to not apply element beam.
			% Apply dipole primary beam correction. 
			% NOTE: Multiply with elembeam pulls up the image edges.
			img.map = img.map .* (elembeam.*imgmask); 
			
			if (wr2file == 0)
				figure (hdl);
	   	     	imagesc(img.l, img.m, (real(img.map).*imgmask));
	   	     	% imagesc((real(img.map))); % Uncomment for pixel axis.
				if isempty (caxisrng) == 0
					caxis (caxisrng);
				end;
				% [sel, sel_l, sel_m] = overplotcat (img.tobs, srclist3CR, ... 
				% 									 10, hdl, true);
	        	set(gca, 'FontSize', 16);
			
		        title(sprintf('Rec: %d, Time:%.2f, Freq:%.2f',offset+ts,...
						img.tobs, img.freq));
		        axis equal
		        axis tight
				% To match orientation with station images
   				set (gca, 'YDir', 'Normal'); 
   				set (gca, 'XDir', 'Reverse'); 
		        ylabel('South \leftarrow m \rightarrow North');
		        xlabel('East \leftarrow l \rightarrow West');
		        set(colorbar, 'FontSize', 16);
			else
				% Write image to output file.
				img.pix2laxis = length (img.l);
				img.pix2maxis = length (img.m);
				wrimg2bin (fimg, img);
			end;
	
			% Read in next visibility set.
			for ind = 1:skip
				[acc, img.tobs, img.freq] = readms2float (fcal, -1, -1, 288);
				rec = readcalsol (fsol);
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

	fclose (fcal); fclose (fsol);
