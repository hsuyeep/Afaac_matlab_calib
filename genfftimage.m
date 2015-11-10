% pep/18Oct12
% Arguments:
%	fname   :   Name of file or matrix containing calibrated or uncalibrated visibilities.
%               NOTE: If fname is a matrix of visibilities, time of observation goes on the
%               end;
%	ntslices:   Number of timeslices to image. -1 => all timeslices. If fname is a matrix,
%               the number of timeslices is the last dimension of the matrix.
%	offset  :   Initial offset in visibility file from which to start imaging.
%               Set to 0 for no offset.
%   skip    :   Number of records to skip, in the interest of faster imaging.
%   integ   :   Number of records over which to integrate (after rephasing
%               to the center of the integration period).
%	posfilename:Name of file containing the local positions of the array 
%               configuration corresponding to the data in fname.
%	mosaic  :   Bool controlling generation of mosaic images. If false,  zenith
%               pointing images (less computation) are generated.
%	caxisrng:   2-element vector with the coloraxis range for display.
%	wr2file :   Variable controlling writing of generated images to file. 
%               0 = Don't write to file.
%               1 = the input filename with 'img' appended is created and written
%                   out using floats. NOTE: If writing to file,images are not shown.
%               2 = A separate fits file per timeslice is created. WCS has a 
%                   a SIN projection.
%   elbeam  :   Bool controlling applying an element beam pattern correction.
%               NOTE: Currently only for X-dipole.
%   pol     :   Polarization to choose: 0 => X, 1 => Y, 2 => I(?)
%  radec    :   Flag to control if RA/DEC plane projected images are to be created.
%  accum    :   Flag controlling the generation of accumulated images.
%               One Accumulated image over the full run.
%  Returns:
%	img_l   :   The l-axis of the generated image.
%	img_m   :   The m-axis of the generated image.
%     img   :   The actual image matrix.
%	NOTE :      The values returned are for the last image in the image 
%				timeseries!


function [img_l, img_m, img, acc_radecmap, acc_localmap] =  ... 
    genfftimage (fname, ntslices, offset, skip, integ, posfilename, mosaic, caxisrng,...
				 wr2file, elbeam, pol, radec, accum, varargin)

    % Default assumption is that fname is a binary file of visibilities
    acc_is_matrix = 0;

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

    % To generate a filename indicating whether PBCORR is applied or not.
	if (elbeam == 0)
		elstr = 'noel';
	else
		elstr = 'el';
	end;

    % For projecting l,m images to RA/DEC space.
	if (radec == 1)
		projimg = figure;
		ra_grid = linspace (0,2*pi, gparm.uvpad)*12/pi; % Convert to hours
		de_grid = linspace (-pi/2,pi/2,gparm.uvpad)*180/pi; % Convert to deg.
		if (accum == 1)
			acc_radecmap = zeros (gparm.uvpad);
			acc_localmap = zeros (gparm.uvpad);
		end;
    else
		acc_radecmap = [];
		acc_localmap = [];
	end;

    % Figure out if the passed fname is a filename of visibilities, or an actual
    % Matrix.
    % Anonymous function to generate char name of variable in workspace.
    % Needed as exist () takes only charstrings.
    namefromvar = @(x) inputname(1);
    if (exist (namefromvar(fname)) == 1) % If fname is a variable in the current workspace.
        acc_is_matrix = 1;
        acm_tseries = fname;
        tobs = varargin{1};
        freq = varargin{2};
    end;

	% Note: we write out the accumulated image even if wr2file == 0
	hdl = figure;
	if (accum == 1)
        if (acc_is_matrix == 0)
		    k = strfind (fname, '.bin');
		    if (mosaic == 1)
	    		accimgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_mosimg_radec.fig');
	    	else
	    		accimgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_fftimg_radec.fig');
	    	end;
        else % Case of fname being a workspace variable
		    if (mosaic == 1)
	    		accimgfname = sprintf ('%d_%d_%s%s', tobs(1), offset, elstr, '_mosimg_radec.fig');
	    	else
	    		accimgfname = sprintf ('%s_%d_%s%s', tobs(1), offset, elstr, '_fftimg_radec.fig');
            end;
        end;

		if (exist (accimgfname, 'file') == 2)
			fprintf (2, 'Overwriting existing file: %s. Continue? (Ctrl-C to kill)\n', accimgfname);
			pause;   % To prevent overwriting already written files!
		end;
	end;

	% Obtain image parameters from visibility file.
    if (acc_is_matrix == 0)
		fin = fopen (fname, 'rb');
		[acc, img.tobs, img.freq] = readms2float (fin, offset, -1, 288);
    else
        acc =  acm_tseries(1,:,:);
        img.tobs = tobs(1);
        img.freq = freq(1);
    end;

	Nelem = size (acc, 1);
	lambda = 299792458/img.freq; 		% in m.
	% duv = lambda/2;

	if (wr2file == 1)
        if (acc_is_matrix == 0)
			% NOTE:  Writing either accumulated image, or individual ones, not both at the same time.
			k = strfind (fname, '.bin');
			if (mosaic == 1)
				imgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_mosimg.bin');
			else
				imgfname = sprintf ('%s_%d_%s%s', fname(1:k-1), offset, elstr, '_fftimg.bin');
			end;
        else
			if (mosaic == 1)
				imgfname = sprintf ('%s_%d_%s%s', tobs(1), offset, elstr, '_mosimg.bin');
			else
				imgfname = sprintf ('%s_%d_%s%s', tobs(1), offset, elstr, '_fftimg.bin');
			end;
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

    elseif (wr2file == 2)
       fimg = -1;
       fits_postfix = sprintf ('_%d_.fits', img.freq); % TODO: TO BE TESTED
	end;

	if (ntslices == -1) 			% Get number of records from the file
        if (acc_is_matrix == 0)
		    [ntslices, tmin, tmax, dt] = getnrecs (fname); 
        else
            ntslices = size (acc,1);
            tmin = tobs(1);
            tmax = tobs(length(tobs));
            dt = tobs(2) - tobs(1);
        end;
	end;

	
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

	% Generate primary beam for given frequency from full EM station level 
	% simulations.
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

	if (isempty (skip)) skip = 0; end;
	fprintf (1, 'Skipping %d recs.\n', skip);
	ts = 1;

    % Figure out if visibility integration is required.
    if (isempty (integ) == 1) integ = 1; end;

    % Generate the datastructure to hold visibilities before integration.
    acm = zeros (size (acc,1), size(acc,2), integ);
    tobs_win = zeros (1, integ);

    % Initialize the integrated visibility window.
    if (acc_is_matrix == 0)
	    for ind = 1:integ
		    [acc, tobs_tmp, img.freq] = readms2float (fin, offset, -1, 288);
	        acm(:,:,ind) = acc;
	        tobs_win(ind) = tobs_tmp;
	    end;
	    ts = ts + integ;
    else
        acm = acm_tseries (:,:,1:integ);
        tobs_win = tobs(1:integ);
    end;

	while (ts < ntslices)
	% for ts = 1:ntslices % Can't accomodate skips within a for loop.
        if (integ > 1)
            [~, acc, img.tobs] = integvis (acm, integ, tobs_win, img.freq, posITRF, 0);
        else
            acc = acm(:,:,1);
            img.tobs = tobs_tmp;
        end;

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
			
			% Accumulated image is shown regardless of wr2file status
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
			elseif (wr2file == 1)
				if (accum == 1)
					saveas (projimg, accimgfname);
				end;
				% Write image to output file.
				img.pix2laxis = length (img.l);
				img.pix2maxis = length (img.m);
				wrimg2bin (fimg, img);
            else
                fitsfname = sprintf ('%d_%s', img.tobs, fits_postfix);
                fimg = fopen (fitsfname, 'wb');
                wrimg2fits (fimg, img);
                fclose (fimg);
			end;
	
			% Skip visibilities
            if (acm_is_matrix == 0)
				for ind = 1:skip
					[acc, img.tobs, img.freq] = readms2float (fin, -1, -1, 288);
				end;
			    ts = ts + skip;

				% Read in next visibility set.
			    for ind = 1:integ
				    [acc, tobs_tmp, img.freq] = readms2float (fin, offset, -1, 288);
			        acm(:,:,ind) = acc;
			        tobs_integ(ind) = tobs_tmp;
			    end;
			    ts = ts + integ;
            else
                acm = acm_tseries(:,:,ts:ts+integ);
			    ts = ts + integ;
            end;
            % 
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

	% Normalize the accumulated images
	if (accum == 1)
		naccum = floor((ts - offset)/skip);
		acc_localmap = acc_localmap/naccum;
		acc_radecmap = acc_radecmap/naccum;
	end;

	% Return out only the last image
	img_l = img.l; img_m = img.m; img = img.map;
    if (acc_is_matrix == 0)
    	fclose (fin);
    end;
	if (wr2file == 1)		
		fclose (fimg);
	end;
