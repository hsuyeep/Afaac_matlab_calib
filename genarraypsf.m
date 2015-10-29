% Script to generate the PSF of a specified AARTFAAC array configuration.
% pep/15Jul13
% Arguments:
%	poslocal: The array element positions in CS002 local coords. 
%   flagant : list of flagged antennas
%	   freq : The frequency of operation, in Hz.
%     tparm : Visibility taper parameters to apply. (See taper.m for fields)
%     wparm : Visibility weighting parameters to apply. (See genvisweight.m for 
%			  fields)
%     gparm : Gridding parameters to apply. Fields as below:
%		 .fft: Bool indicating FFT or DFT imaging.
%        .duv: The grid spacing in meters.
%        .Nuv: The grid size in pixels.
%        .uvpad : Padding to apply for higher spatial resolution in 2D FFT.
%       deb : Bool turning on debug information.

% Returns:
% l,m    : local coordinates for this PSF as a matrix.
%	psf  : Actual PSF as an image matrix
%  weight: The weighing matrix used.
%  intap : The applied inner taper for this psf.
%  outtap: The applied outer taper for this psf.
% uvdist : The uvdistance of each visibility from the UV plane center.
%    vis : The gridded and weighted visibility.

function [l, m, psf, weight, intap, outtap,uvdist, vispad] = ...
	genarraypsf (posfilename, flagant, freq,  tparm, wparm, gparm, deb)
	% Default values are for 60MHz, and image just the full Fov (-1<l<1). 
	if (isempty (gparm) == 1)
		gparm.type = 'pillbox';
		% gparm.type = 'gaussian';
		gparm.duv = 0.5;  % Now in wavelengths
		gparm.Nuv = 241;
		gparm.lim = 1.5; % In wavelength
		gparm.uvpad = 256; 
		gparm.fft = 0;
		gparm.pa(1) = 0.5; gparm.pa(2) = gparm.pa(1);
	end;

	% Load coordinates
	load (posfilename, 'poslocal', 'posITRF');

	% Flag antennas as specified in flagant;
	if (~isempty (flagant))
		sel = zeros (length (posITRF), 1);
		sel (flagant) = 1;
		posITRF = posITRF(sel == 0,:);
		poslocal = poslocal (sel == 0, :);
	end;

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
		tparm.pa = [-1 -1 -1 -1]; weighting_noise = 0;
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

	% Create the uloc, vloc meshgrid for frequency independent coords.
	uloc_flag = meshgrid (poslocal (:,1)) - meshgrid (poslocal(:,1)).';
	vloc_flag = meshgrid (poslocal (:,2)) - meshgrid (poslocal(:,2)).';
	
%%%% Now flagging earlier!
%	if (isempty(flagant) == 0)
%		[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant);
%		fprintf (2, 'Flagging %d antennas.\n', length (flagant));
%	else 
%		uloc_flag = uloc; vloc_flag = vloc;
%	end;

	% Set the visibilities at corresponding (uv) points to 1. add their weights
	% and taper.
	weight (weight == 0) = 1;
	acc = ones (length (weight), 1) .* 1./weight .* intap .* outtap;

	% Rephase the ACM to the required pointing
	% [re_acm, rot_uvw] = rephaselm (acc, lreq, mreq, tobs, freq, poslocal);
	if (gparm.fft == 1)
		fprintf (1, '-- Generating FFT image, gridding with supplied params.\n');
		[rdsky, psf, vispad, l, m] = fft_imager_sjw_radec (acc, uloc_flag, ...
						vloc_flag, gparm, [], [], 0, freq, 0);
	else
		fprintf (1, '-- Generating DFT image.\n');
		l = linspace (-1, 1, gparm.uvpad); m = l;
		% l = [-1:0.01:1]; m = l;
		nelem = sqrt (length(acc));
		psf = acm2skyimage (reshape (acc, [nelem nelem]), poslocal(:,1), poslocal(:,2), freq, l, m);
		vispad = fftshift (fft2(psf));
	end;

	% Generate weighting noise for this configuration
	% see SIRA-2, pp. 131 for the formula
	weighting_noise = sqrt(sum((intap .* outtap).^2.*(1./weight).^2))./sum((intap.*outtap).*(1./weight)); 
	fprintf (2, 'Weighting noise for this arrangement: %f.\n', weighting_noise);

	% If required, generate plots of PSF.
	if (deb > 0)
		subplot (122);
		d = intersect (find (l > -0.2), find (l < 0.2)); 
		% mesh (l(d), m(d), double (20*log10 (abs (psf(d, d))/max(max(abs(psf)))))); colorbar;
		% xlabel ('m'); ylabel ('l'); zlabel ('dB');
		% axis ([-1 1 -1 1 -60 0]);
		imagesc (l(d), m(d), 20*log10 (abs (psf(d,d))/max(max(abs(psf))))); colorbar;
		caxis ([-75 0]); % Match colorbar limit to line plot limit.
		% title (sprintf ('AARTFAAC PSF: weight=%s, cellrad=%d ', wparm.type, wparm.cellrad));
		% title (sprintf ('AARTFAAC PSF:')); 

		subplot (221);
		scan = abs (psf (gparm.uvpad/2, d));
		plot (l(d), 20*log10 (scan/max(scan)));
		ylim ([-75 0]);
		% title (sprintf ('l-scan at m=0, taper=%s, weight=%s', taper, weight));
		grid on;
		xlabel ('l'); ylabel ('Power (dB)');

		subplot (223);
		scan = abs (psf (d, gparm.uvpad/2));
		plot (m(d), 20*log10 (scan/max(scan)));
		ylim ([-75 0]);
		% axis ([-0.2 0.2 -60 0]);
		grid on;
		% title (sprintf ('m-scan at l=0, taper=%s, weight=%s', taper, weight));
		xlabel ('m'); ylabel ('Power (dB)');
		% print (gcf, 'psf.eps', '-depsc', '-r300'); 
		mtit (sprintf ('out: %.2f, in: %.2f, wnoise: %.4g', tparm.pa(3), tparm.pa(1), weighting_noise), 'yoff', 0.025);

		figure;
		subplot (211);
		plot (uvdist, acc, '.');
		grid on; axis tight;
		xlabel ('uvdist (m)'); ylabel ('Final weights');

		subplot (212);
		paduvlim = floor (gparm.uvpad/2)*gparm.duv;
    	padgridrng = linspace (-paduvlim, paduvlim, gparm.uvpad);
		ugrid = meshgrid (padgridrng);
		paduvcoor = [ugrid(:) ugrid(:)];
		griduvdist = sqrt (sum(paduvcoor.^2, 2));
		if (gparm.fft == 1)
			[n,x]=hist (griduvdist.*vispad(:), 50);
			% bar (x(2:end), n(2:end));
			plot(griduvdist, vispad(:), '.');
			xlabel ('uvdist (\lambda)');
			ylabel ('visibility weights');
		end;
		
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

