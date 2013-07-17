% Script to generate the PSF of a specified AARTFAAC array configuration.
% pep/15Jul13
% Arguments:
%	poslocal: The array element positions in CS002 local coords. 
%   flagant : list of flagged antennas
%	   freq : The frequency of operation, in Hz.
%     taper : Visibility taper to apply.
%     weight: Visibility weighting to apply.
%        duv: The grid spacing in meters.
%        Nuv: The grid size in pixels.
%     debug : Bool turning on debug information.

% Returns:
%	psf  : Actual PSF as an image matrix

function [psf] = genarraypsf (poslocal, flagant, freq,  taper, weight, duv,...
							  Nuv, debug)
% function [psf] = genarraypsf (poslocal, flagant, freq,  taper, weight, duv,...
% 							  Nuv, debug)

    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 1000;						% size of gridded visibility matrix
    uvpad = 1024;					% specifies if any padding needs to be added

	% Generate the visibility taper function
%	switch lower (taper)
%		case {'gaussian'}
%			fprintf (2, 'Taper mode Gaussian not implemented!\n');
%
%		case {'triangle'}
%			fprintf (2, 'Taper mode Triangle not implemented!\n');
%
%		case {'hanning'}
%			fprintf (2, 'Taper mode Hanning not implemented!\n');
%
%		case {'hamming'}
%			fprintf (2, 'Taper mode Hamming not implemented!\n');
%
%		case {'blackman'}
%			fprintf (2, 'Taper mode Blackman not implemented!\n');
%
%		otherwise,
%			fprintf (2, 'Taper mode %s not known! Not applying taper');
%	end;
%
%	% Generate the visibility density weighting function
%	switch weight
%		case {'uniform'}
%			fprintf (2, 'Weighting mode Uniform not implemented!\n');
%
%		case {'natural'}
%			fprintf (2, 'Weighting mode Natural not implemented!\n');
%
%		otherwise,
%			fprintf (2, 'Taper mode %s not known! Not applying taper');
%	end;
%
%	% Choose gridding interpolation method
%	switch gcf
%		case {'linear'}
%			fprintf (1, 'Applying linear interpolation\n');
%		otherwise,
%			fprintf (2, 'GCF %s not known! Applying linear interpolation..\n');
%	end;
		
	% Create the uloc, vloc meshgrid for frequency independent coords.
	load (posfilename, 'poslocal');
	uloc = meshgrid (poslocal (:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal (:,2)) - meshgrid (poslocal(:,2)).';
	
	[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
	% Set the visibilities at corresponding (uv) points to 1.
	acc = ones (size (uloc_flag));

	[rdsky, psf, vispad, l, m] = fft_imager_sjw_radec (acc, uloc_flag, vloc_flag, duv, Nuv, uvpad, 0, freq, 0);

	% Generate gridded visibilities
	% gridacc = visgrid (acc, duv, uloc_flag, vloc_flag, Nuv);

	% Add the Taper and weighting functions to this dataset, if specified.
	% If not, simply use the uniform weighting and no taper.
	% TODO.

	% Cary out the FFT of the sampled UV place.
	% psf = fft2 (gridacc);

	% If required, generate plots of PSF.
	if (debug > 0)
		subplot (122);
		% imagesc (l, m, 20*log10 (abs (psf))); colorbar;
		mesh (l, m, double (20*log10 (abs (psf)/max(max(abs(psf)))))); colorbar;
		title ('AARTFAAC PSF: ');
		xlabel ('m'); ylabel ('l');
		axis ([-1 1 -1 1 -60 0]);

		subplot (221);
		scan = abs (psf (Nuv/2, :));
		plot (l, 20*log10 (scan/max(scan)));
		axis ([-1 1 -60 0]);
		% title (sprintf ('l-scan at m=0, taper=%s, weight=%s', taper, weight));
		xlabel ('l'); ylabel ('Power (dB)');

		subplot (223);
		scan = abs (psf (:, Nuv/2));
		plot (m, 20*log10 (scan/max(scan)));
		axis ([-1 1 -60 0]);
		% title (sprintf ('m-scan at l=0, taper=%s, weight=%s', taper, weight));
		xlabel ('m'); ylabel ('Power (dB)');
	end; 

