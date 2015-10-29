% Script to carry out pixel level correlations between generated images.
% The resulting correlation function gives information on :
%	- The independence of time sampled images, useful for constraining 
%	  transient rates.
%	- The confusion noise limit, both for classical and Sidelobe confusion.
% Arguments:
%	img      : Image timeseries binary filename, as generated by gen
%	start_rec: Record offset in file to beginning of timeseries.
%	skip     : Number of records to skip in generating the time series.
%	nrecs	 : Number of records in analysis timeseries
%	thresh   : Sigma cutoff for noise estimation.
%	region	 : [center, radius], region of the image on which to operate, in l,m units.
%      deb   : Debug flag, set to 1 for intermediate results.

function [acf_pix] = genpixcorr (fimg, start_rec, skip, nrecs, thresh, region, deb)
	% Check if parameters are valid
	assert (isempty (fimg) == 0);
	if (isempty (start_rec)) start_rec = 1; end;
	if (isempty (skip)) skip = 1; end;
	if (isempty (thresh)) thresh = 3; end; % Default sigma for cutoff.

	% Display information from file.
	fprintf (1, '<-- Operating on file %s.\n', fimg);
	fid = fopen (fimg, 'r');
	img = readimg2bin (fid, 0); 
	fprintf (1, '<-- Image timeseries start : %s, freq: %f.\n', datestr(mjdsec2datenum (img.tobs)), img.freq);
	fprintf (1, '<-- %04dx%04d image, l range : [%.3f, %.3f], m range : [%.3f, %.3f]\n', img.pix2laxis, img.pix2maxis, min(img.l), max(img.l), min(img.m), max(img.m));

	fprintf (1, '<-- Generating desired mask...\n');
	mask = ones (img.pix2laxis, img.pix2maxis);

	[l, m] = meshgrid (img.l, img.m);
	% Generate mask for the specified region, ignoring offset of the center for now.
	mask (l.^2 + m.^2 > region.radius) = 0;

	% Create data structure for storing timeseries of valid pixels.
    fprintf (1, '<-- Pixels in selected region: %d\n', sum(mask(:)));
	pixtseries = zeros (nrecs, sum(mask(:)));
    timg =zeros (nrecs); % To store image times.


	for imgind = 1:nrecs
        fprintf (1, '<-- Collecting image %04d\n', imgind);
		
		% Apply sigma clipping on the masked region, eliminate bright pixels.
		valpix = img.map (mask == 1);
		[mu, v, sel] = robustmean (valpix, thresh);

		% Set the unused pixels to NaNs, so that a consistent matrix can be used.
		valpix (sel == 0) = 0; 

		% store valid pixels
		pixtseries (imgind, :) = valpix;
        timg (imgind) = img.tobs;

        % Read in the next image
    	img = readimg2bin (fid, 0); 
	end;

    if (deb == 1)
        % Show the contents of the valid pixels
        tmpimg = zeros (size(img.map));
        for ind = 1:nrecs
            tmpimg (mask == 1) = pixtseries(ind, :);
            subplot (211);
            imagesc (img.l, img.m, tmpimg); colorbar();
            title (sprintf ('Time: %s', datestr(mjdsec2datenum(timg(ind)))));
            subplot (212);
            hist (pixtseries(ind,:), 50);
            pause();
        end;
    end;

	% 1. Temporal correlation of every pixel: Generate the Auto Correlation Function
	%	 of the timeseries of each pixel. Note, the timeseries is stored in a column.
    fprintf (1, '<-- Calculating temporal correlation...\n');
	Pix = fft(pixtseries,2*nrecs-1); % FFT carried out independently per column
	acf_pix = fftshift (ifft(Pix.*conj(Pix)),1);
	acf_pix = acf_pix ./ repmat (max(acf_pix), 2*nrecs-1, 1); % Normalize corr. coef. to 1.
