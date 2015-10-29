% Function to generate an image from calibrated visibilities.
% pep/18Jul12
% Arguments: 
%     acc   : Complex, Hermitean symmetric Array Correlation Matrix
%     u/v   : u/v coordinates of antenna elements in ITRF coordinates (meters) 
%             wrt. CS002
%    
%  gparm    : -- Gridding parameters, passed onto genvisgrid.m --
%    .type  : Type of grid convolution function.
%    .duv   : Grid size in meters, in the uv domain. Keep at least a couple of 
%             gridpoints per wavelength.
%    .Nuv   : Number of gridpoints sampling the available UV plane.
%    .uvpad : Padded size of uv grid, for imaging with higher resolution than
%             warranted by uv coverage.
%    .lim   : Limit (in meters) of the GCF. Use if appropriate to GCF type.
%	 .pa(1:4): Additional parameters for describing	the GCF.

%  wparm    : -- Weighting parameters, passed onto genvisweight.m -- 

%  tparm    : -- Tapering parameters, passed onto taper.m -- 
% 
%     t_obs : Time of observation in MJD seconds, for converting local 
%             coordinates to absolute RA/DEC coordinates.
%     freq  : Frequency of observation in Hz, of this ACM.
%     radec : Switch to turn on (=1) or off(=0) the generation of images 
%             reprojected onto the RA/DEC plane.
% Returns: 
%    radecskymap: Reprojected skymap in RA/DEC coordinates, 
%    lmskymap   : Map in lm local coordinates.
%    vispad     : Gridded and padded visibilities
%    l,m        : l,m coordinates for this skymap.

function [radecskymap, lmskymap, vispad, l, m] =  ... 
		fft_imager_sjw_radec(acc, u, v, gparm, wparm, tparm, t_obs, freq, radec)

%  Example gridding parameters
%	gparm.type = 'wijnholds';
%	gparm.duv = duv;
%	gparm.Nuv = Nuv;
%	gparm.uvpad = uvsize;
%	gparm.lim = 0;
%	gparm.pa = [0 0 0 0];
	[vispad, gridviscnt, padgridrng, kern] = genvisgrid (acc, u, v, gparm, freq, 0);
    
    % compute image
    % ac = zeros(size(vispad));
    % ac(256, 256) = 100;
    % vispad = vispad + ac;
    
    % FFT imaging: Just take FFT of the gridded visibilities.
    vispad = conj(flipud(fliplr(fftshift(vispad))));
    skymap = fftshift(fft2(vispad));
    
    % Create l,m axis corresponding to choices of duv
    % NOTE: Resolution of image is determined by total array aperture extent, 
    % resolution = lambda/D. 
    % dl = (299792458/(freq * gparm.uvpad * gparm.duv)); % dimensionless, in dir. cos. units
    dl = (1/(gparm.uvpad * gparm.duv)); % dimensionless, in dir. cos. units
	% Changed due to gparm.duv now being specified in wavelengths, not meters.
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * gparm.uvpad/ 2;
	l = linspace (-lmax, lmax, gparm.uvpad);
    % l = [-lmax:dl:lmax-1e-3];
    m = l;  % Identical resolution and extent along m-axis
    
    % Create a mask to mask out pixels beyond the unit circle (these are below 
	% the horizon.)
    mask = NaN (length(l));
    mask(meshgrid(l).^2 + meshgrid(l).'.^2 < 1) = 1;
    lmskymap = single (real(skymap) .* mask)';
	% Transpose added on 03Mar15, such that images with axis match the real sky.
	% Rescaling required due to the unnormalized FFT carried out by fft2.
    % lmskymap = single ((skymap) .* mask)'./sqrt(length(acc));
    % disp (['-->Max/min from lm skymap: ' num2str(max(max(lmskymap))) ' ' ... 
    % 		num2str(min(min(lmskymap)))]);
    
    % --- Section converting image from local coordinates to RA/DEC coordinates
    if (radec == 1)
      disp ('Converting to RA/DEC image');

	  % Holding area for pixels converted from lm to RA/Dec. coords.
	  % Regular sampling in lm coordinates leads to irregular sampling in RA/Dec
      alpha = zeros (gparm.uvpad, 1);
      delta = zeros (gparm.uvpad, 1);
      % t_obs  = t_obs/86400 + 2400000.5; % Convert to MJD day units

	  % Convert image coordinates to RA/Dec. coordinates.
      [alpha, delta] = lmtoradec (l, l, JulianDay (mjdsec2datenum (t_obs)), 6.869837540, 52.915122495);
      sel = ~isnan (alpha(:));
      
	  % Create a Matlab interpolation object,
      radecimage = TriScatteredInterp (alpha (sel), delta (sel), ... 
										real(skymap (sel)));

      % Create the regularly sampled RA/dec plane
      [ragrid, decgrid] = meshgrid (linspace (0,2*pi, gparm.uvpad), ...
									linspace (-pi/2,pi/2, gparm.uvpad));
          
      % generate samples from the interpolated model skyimage
      radecskymap = radecimage (ragrid, decgrid);
    else 
      % disp ('Not converting to RA/DEC images');
      radecskymap = lmskymap;
    end
