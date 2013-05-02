% Script to generate FFT images from calibrated visibilities available in a 
% .bin file.
% pep/18Oct12
% Arguments:
%	fname : 	Name of file containing calibrated or uncalibrated visibilities.
%	ntslices:   Number of timeslices to image. -1 => all timeslices.
%	offset:		Initial offset in visibility file from which to start imaging.
%				Set to 0 for no offset.
%	posfilename:Name of file containing the local positions of the array 
%				configuration corresponding to the data in fname.
%	mosaic:		Bool controlling generation of mosaic images. If false,  zenith
%				pointing images (less computation) are generated.
%	caxisrng: 	2-element vector with the coloraxis range for display.
%	wr2file:Bool controlling writing of generated images to file. If true,
%				the input filename with 'img' appended is created and written
%				out using floats. NOTE: If writing to file, images are not shown.

function [img_l, img_m, img] =  ... 
    genfftimage (fname,ntslices, offset, posfilename, mosaic, caxisrng, wr2file)
    % genfftimage (fname,ntslices, offset, posfilename, weight, uvcellsize, mosaic, caxisrng, wr2file)
	radec = 0;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500; %1000                % size of gridded visibility matrix
    uvpad = 512; %1024              % specifies if any padding needs to be added
	nfacet = 3;
	facetsize = 256;

	if (wr2file == 0)
		hdl = figure;
	else
		k = strfind (fname, '.bin');
		if (mosaic == 1)
			imgfname = [fname(1:k-1) '_mosimg.bin'];
		else
			imgfname = [fname(1:k-1) '_fftimg.bin'];
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
	duv = lambda/2;
	
	station = 0;
	% For imaging from a single station
	% acc_tmp = acc;
	% acc = zeros (48, 48);
	% acc = acc_tmp (1:48, 1:48);
	% acc = acc_tmp (49:96, 49:96);
	% acc = acc_tmp (97:144, 97:144);
    dl = (299792458/(img.freq * uvpad * duv)); % dimensionless, in dir. cos. units
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * uvpad / 2;
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
					duv, Nuv, uvpad, img.tobs, img.freq, radec);
	else
		% Image current timeslice. Generate a mosaic image.
   		[img.map, img.l, img.m] = ... 
				genmosaic(acc, img.tobs, img.freq, nfacet, facetsize, ...
						uloc(:), vloc(:), posITRF, 0);
	end;
	mask = zeros (size (img.map));
	mask (meshgrid (img.l).^2 + meshgrid(img.m).'.^2 < 0.9) = 1;
	disp ('NOTE! NOTE! Working only with LBA X-dipoles primary beam correction now!');
	addpath 'LBA_beam/CS1/';
	% elembeam = calculateLBAbeam (img.l, img.m, img.freq, [1:2:96]); 
	elembeam = ones (size (mask));
	elembeam = max(max(elembeam(:,:,1))) ./ elembeam (:,:,1);

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

	for ts = 1:ntslices

		if (isempty (acc) == false)
			fprintf (1, 'Slice: %3d, Time:%.2f, Freq:%.2f.\n', ... 
					 offset+ts,img.tobs, img.freq); 
			acc_weighted = acc .* weightmask;

			if (mosaic == 0)
				% Image current timeslice. Generate a zenith image.
		   		[radecmap, img.map, calvis, img.l, img.m] = ... 
				  fft_imager_sjw_radec (acc_weighted(:), uloc(:), vloc(:), ... 
							duv, Nuv, uvpad, img.tobs, img.freq, radec);
			else
				% Image current timeslice. Generate a mosaic image.
		   		[img.map, img.l, img.m] = ... 
						genmosaic(acc_weighted, img.tobs, img.freq, nfacet, facetsize, ...
								uloc(:), vloc(:), posITRF, 0);
			end;

			% Apply dipole primary beam correction. 
			% Comment to not apply element beam.
			img.map = img.map .* (elembeam.*mask); 
			
			if (wr2file == 0)
				figure (hdl);
	   	     	imagesc(img.l, img.m, (real(img.map)));
				if isempty (caxisrng) == 0
					caxis (caxisrng);
				end;
				% [sel, sel_l, sel_m] = overplot3cr (img.tobs, srclist3CR, ... 
				%									 100, hdl);
	        	set(gca, 'FontSize', 16);
			
		        title(sprintf('Rec: %d, Time:%.2f, Freq:%f',offset+ts,img.tobs, ...
						img.freq));
		        axis equal
		        axis tight
		        xlabel('South \leftarrow m \rightarrow North');
		        ylabel('East \leftarrow l \rightarrow West');
		        set(colorbar, 'FontSize', 16);
			else
				% Write image to output file.
				img.pix2laxis = length (img.l);
				img.pix2maxis = length (img.m);
				wrimg2bin (fimg, img);
			end;
	
			% Read in next visibility set.
			[acc, img.tobs, img.freq] = readms2float (fin, -1, -1, 288);
			if (isempty (acc))
				disp ('End of file reached!');
				break;
			end;
			% pause (0.01);
			% pause;
		else
			disp ('File end reached!');
		end;
	end;
	fclose (fin);
	if (wr2file == 1)		
		fclose (fimg);
	end;
