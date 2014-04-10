% Script to generate the PSF of a specified AARTFAAC array configuration.
% pep/15Jul13
% Arguments:
%	poslocal: The array element positions in CS002 local coords. 
%   flagant : list of flagged antennas
%	   freq : The frequency of operation, in Hz.
%     tparm : Visibility taper parameters to apply. (See taper.m for fields)
%     wparm : Visibility weighting parameters to apply. (See genvisweight.m for 
%			  fields)
%     iparm : Imaging parameters to apply. Fields as below:
%		 .fft: Bool indicating FFT or DFT imaging.
%        .duv: The grid spacing in meters.
%        .Nuv: The grid size in pixels.
%        .uvpad : Padding to apply for higher spatial resolution in 2D FFT.
%       deb : Bool turning on debug information.

% Returns:
% l,m    : local coordinates for this PSF as a matrix.
%	psf  : Actual PSF as an image matrix
%  weight: The weighing matrix used.
%  taper : The applied taper for this psf.

function [l, m, psf, weight, intap, outtap] = ...
	genarraypsf (posfilename, flagant, freq,  tparm, wparm, iparm, deb)
	% Default values are for 60MHz, and image just the full Fov (-1<l<1). 
	if (isempty (iparm) == 1)
		iparm.duv = 2.5; 
		iparm.Nuv = 1000;
		iparm.uvpad = 1024; 
		iparm.fft = 1;
	end;

	% Load coordinates
	load (posfilename, 'poslocal', 'posITRF');

	% Generate an appropriate taper
	if (isstruct (tparm))
		[intap, outtap, den, mask, uvdist] =  taper (posITRF, tparm, -1, freq, deb);
	else
%		tparm.type = 'Gaussian';
%		tparm.minlambda = 10;    % Spatial filter lower cutoff, Units of lambda. 
%		tparm.maxmeters = 350;	 % Spatial filter upper cutoff, Units of meters.
%		tparm.pa(1) = -1; %0.2;	 % Inner gaussian taper sigx, Units of lambda. 
%								 % -1 => no inner taper. 
%		tparm.pa(2) =tparm.pa(1);% Inner gaussian taper sigy,Units of lambda
%		tparm.pa(3) = 100; 	     % Outer gaussian taper sigx, Units of lambda
%		tparm.pa(4) =tparm.pa(3);% Outer gaussian taper sigy, Units of lambda
%		[intap, outtap, den, mask, uvdist] =  taper (posITRF, tparm, -1, freq, deb);
		
		fprintf (2, '!! taper unspecified! Continuing with no tapering.\n');
		intap = ones (size (posITRF,1)*size(posITRF,1), 1); outtap = intap;
		u = meshgrid(posITRF(:, 1)) - meshgrid(posITRF(:, 1)).';
	   	v = meshgrid(posITRF(:, 2)) - meshgrid(posITRF(:, 2)).';
	    w = meshgrid(posITRF(:, 3)) - meshgrid(posITRF(:, 3)).';
	    uvw = [u(:), v(:), w(:)];
        normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	    uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
	end;

	% Generate the visibility density weighting function
	if (isstruct (wparm))
		weight = genvisweight (posITRF, wparm, 1);
	elseif isempty(wparm)
		fprintf (2, '!! wparm empty! Continuing with Natural weighting.\n');
		weight = ones (length(uvdist), 1);
	else
		fprintf (2, '!! wparm unspecified! Continuing with user specified weighting.\n');
		weight = wparm;
	end;

%
%	% Choose gridding interpolation method
%	switch gcf
%		case {'linear'}
%			fprintf (1, 'Applying linear interpolation\n');
%		otherwise,
%			fprintf (2, 'GCF %s not known! Applying linear interpolation..\n');
%	end;
		
	% Create the uloc, vloc meshgrid for frequency independent coords.
	uloc = meshgrid (poslocal (:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal (:,2)) - meshgrid (poslocal(:,2)).';
	
	if (isempty(flagant) == 0)
		[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
		fprintf (2, 'Flagging %d antennas.\n', length (flagant));
	else 
		uloc_flag = uloc; vloc_flag = vloc;
	end;

	% Set the visibilities at corresponding (uv) points to 1. add their weights
	% and taper.
	weight (weight == 0) = 1;
	acc = ones (length (weight), 1) .* 1./weight .* intap .* outtap;

	if (iparm.fft == 1)
		[rdsky, psf, vispad, l, m] = fft_imager_sjw_radec (acc, uloc_flag, ...
						vloc_flag, iparm.duv, iparm.Nuv, iparm.uvpad, 0, freq, 0);
	else
		l = linspace (-1, 1, iparm.Nuv); m = l;
		% l = [-1:0.01:1]; m = l;
		psf = acm2skyimage (reshape (acc, [288 288]), poslocal(:,1), poslocal(:,2), freq, l, m);
	end;

	% Generate gridded visibilities
	% gridacc = visgrid (acc, duv, uloc_flag, vloc_flag, Nuv);


	% Generate weighting noise for this configuration
	% see SIRA-2, pp. 131 for the formula
	weighting_noise = sqrt(sum((intap .* outtap).^2.*(1./weight).^2))./sum((intap.*outtap).*(1./weight)); 
	fprintf (2, 'Weighting noise for this arrangement: %f.\n', weighting_noise);

	% If required, generate plots of PSF.
	if (deb > 0)
		subplot (122);
		% imagesc (l, m, 20*log10 (abs (psf))); colorbar;
		mesh (l, m, double (20*log10 (abs (psf)/max(max(abs(psf)))))); colorbar;
		% title (sprintf ('AARTFAAC PSF: weight=%s, cellrad=%d ', wparm.type, wparm.cellrad));
		title (sprintf ('AARTFAAC PSF:')); 
		xlabel ('m'); ylabel ('l'); zlabel ('dB');
		axis ([-1 1 -1 1 -60 0]);

		subplot (221);
		scan = abs (psf (iparm.Nuv/2, :));
		plot (l, 20*log10 (scan/max(scan)));
		axis ([-0.2 0.2 -60 0]);
		% axis ([-1 1 -60 0]);
		% title (sprintf ('l-scan at m=0, taper=%s, weight=%s', taper, weight));
		grid on;
		xlabel ('l'); ylabel ('Power (dB)');

		subplot (223);
		scan = abs (psf (:, iparm.Nuv/2));
		plot (m, 20*log10 (scan/max(scan)));
		axis ([-0.2 0.2 -60 0]);
		grid on;
		% title (sprintf ('m-scan at l=0, taper=%s, weight=%s', taper, weight));
		xlabel ('m'); ylabel ('Power (dB)');
		% print (gcf, 'psf.eps', '-depsc', '-r300'); 

		figure;
		plot (uvdist, acc, '.');
		grid on; axis tight;
		xlabel ('uvdist (m)'); ylabel ('Final weights');
		% print (gcf, 'weights.eps', '-depsc', '-r300'); 
		% title (sprintf ('Weight: %s, cellrad: %d', wparm.type, wparm.cellrad));
%		tparm.type = 'Gaussian';
%		tparm.minlambda = 10;    % Spatial filter lower cutoff, Units of lambda. 
%		tparm.maxmeters = 350;	 % Spatial filter upper cutoff, Units of meters.
%		tparm.pa(1) = -1; %0.2;	 % Inner gaussian taper sigx, Units of lambda. 
%								 % -1 => no inner taper. 
%		tparm.pa(2) =tparm.pa(1);% Inner gaussian taper sigy,Units of lambda
%		tparm.pa(3) = 100; 	     % Outer gaussian taper sigx, Units of lambda
%		tparm.pa(4) =tparm.pa(3);% Outer gaussian taper sigy, Units of lambda
	end; 

