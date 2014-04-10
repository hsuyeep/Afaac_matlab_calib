% Function to generate an image from calibrated visibilities.
% pep/18Jul12
% Arguments: 
%     acc   : Complex, Hermitean symmetric Array Correlation Matrix
%     u/v   : u/v coordinates of antenna elements in ITRF coordinates (meters) 
%             wrt. CS002
%     duv   : Grid size in meters, in the uv domain. Keep at least a couple of 
%             gridpoints per wavelength.
%     Nuv   : Number of gridpoints sampling the available UV plane.
%     uvsize: Padded size of uv grid, for imaging with higher resolution than
%             warranted by uv coverage.
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
		fft_imager_sjw_radec(acc, u, v, duv, Nuv, uvsize, t_obs, freq, radec)

%
%    % create object for interpolation
%    vis = zeros(Nuv);
%	missed_vis = 0;   % cumulative count of ignored visibilities due to 
%					  % gridded value exeeding grid size.
%    %W = zeros(Nuv);
%    for idx = 1:length(u(:))            % For every recorded visibility
%    
%    	% Get amp. and direction vector of visibility
%        ampl = abs(acc(idx));			
%        phasor = acc(idx) / ampl;
%    
%    	% Determine the grid along U-axis in which observed visibility falls.
%        uidx = u(idx) / duv + Nuv / 2;  
%        uidxl = floor(uidx);		% Find the lower and higher gridded U-value.
%        uidxh = ceil(uidx);
%    
%        % Find absolute distance of measured visibility from grid points, in the
%    	% U-direction only.
%        dul = abs(uidx - uidxl);		
%        duh = abs(uidx - uidxh);
%    
%        % Distribute the visiblity amplitude among the two grid points on the 
%		% U-axis in proportion to their distance from the observed visiblity.
%        sul = duh * ampl;
%        suh = dul * ampl;
%        
%    	% Determine the grid along V-axis in which observed visibility falls.
%        vidx = v(idx) / duv + Nuv / 2;
%        vidxl = floor(vidx);		% Find the lower and higher gridded V-value.
%        vidxh = ceil(vidx);
%    
%        % Find absolute distance of measured visibility from grid points, in the
%    	% V-direction only.
%        dvl = abs(vidx - vidxl);
%        dvh = abs(vidx - vidxh);
%    
%    	% Distribute the HIGHER u-grid point's share of the observed visibility
%		% amp. between the higher and lower V-grid point.
%        sull = dvh * sul;
%        suhl = dvh * suh;
%    
%    	% Distribute the LOWER u-grid point's share of the observed visibility 
%		% amp. between the higher and lower V-grid point.
%        sulh = dvl * sul;
%        suhh = dvl * suh;
%        
%    	% Now that the observed visiblity amplitude is distributed among its 
%    	% surrounding 4 grid points, fill the gridded visibility matrix with
%    	% vectors with the same phase as the original observed visibility.
%    	% NOTE: Adding the 4 vectors at the corners of the grid square will give
%    	% back the original ungridded observed visibility.
%    	% NOTE: We need to accumulate this to prevent overwriting the gridded 
%		% values from a visibility from a neighbouring grid square.
%		% fprintf (1, 'bline: %06d (%03.2f,%6.2f), uidx: %03d %03d, vidx: %03d, %03d\n', idx, u(idx), v(idx), uidxl, uidxh, vidxl, vidxh);
%		if ((uidxl < 1) || (uidxh < 1) || (uidxl > Nuv) || (uidxh > Nuv))
%			missed_vis = missed_vis + 1;
%			continue;
%		end;
%
%		if ((vidxl < 1) || (vidxh < 1) || (vidxl > Nuv) || (vidxh > Nuv))
%			missed_vis = missed_vis + 1;
%			continue;
%		end;
%
%        vis(uidxl, vidxl) = vis(uidxl, vidxl) + sull * phasor;
%        vis(uidxl, vidxh) = vis(uidxl, vidxh) + sulh * phasor;
%        vis(uidxh, vidxl) = vis(uidxh, vidxl) + suhl * phasor;
%        vis(uidxh, vidxh) = vis(uidxh, vidxh) + suhh * phasor;
%        
%        %W(uidx, vidx) = W(uidx, vidx) + 1;
%    end
%
%	if (missed_vis > 0)
%	 	fprintf (2, 'Missed vis: %d\n', missed_vis); 
%	end;
%
%    % zero padding to desired (u,v)-size
%    N = size(vis, 1);
%    N1 = floor((uvsize - N) / 2);
%    N2 = ceil((uvsize + 1 - N) / 2) - 1;
%    
%    % Surround gridded visibilities with 0-padding to create padded visibility 
%    % matrix.
%    vispad = [zeros(N1, uvsize); ...
%              zeros(N, N1), vis, zeros(N, N2); ...
%              zeros(N2, uvsize)];
%    vispad(~isfinite(vispad)) = 0;

	parm.type = 'wijnholds';
	parm.duv = duv;
	parm.Nuv = Nuv;
	parm.uvpad = uvsize;
	parm.lim = 0;
	parm.pa = [0 0 0 0];
	vispad = genvisgrid (acc, u, v, parm, 0);
    
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
    dl = (299792458/(freq * uvsize * duv)); % dimensionless, in dir. cos. units
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * uvsize / 2;
	l = linspace (-lmax, lmax, uvsize);
    % l = [-lmax:dl:lmax-1e-3];
    m = l;  % Identical resolution and extent along m-axis
    
    % Create a mask to mask out pixels beyond the unit circle (these are below 
	% the horizon.)
    mask = NaN (length(l));
    mask(meshgrid(l).^2 + meshgrid(l).'.^2 < 1) = 1;
    % lmskymap = single (real(skymap) .* mask);
    lmskymap = single ((skymap) .* mask);
    % disp (['-->Max/min from lm skymap: ' num2str(max(max(lmskymap))) ' ' ... 
    % 		num2str(min(min(lmskymap)))]);
    
    % --- Section converting image from local coordinates to RA/DEC coordinates
    if (radec == 1)
      disp ('Converting to RA/DEC image');

	  % Holding area for pixels converted from lm to RA/Dec. coords.
	  % Regular sampling in lm coordinates leads to irregular sampling in RA/Dec
      alpha = zeros (uvsize, 1);
      delta = zeros (uvsize, 1);
      % t_obs  = t_obs/86400 + 2400000.5; % Convert to MJD day units

	  % Convert image coordinates to RA/Dec. coordinates.
      [alpha, delta] = lmtoradec (l, l, JulianDay (mjdsec2datenum (t_obs)), 6.869837540, 52.915122495);
      sel = ~isnan (alpha(:));
      
	  % Create a Matlab interpolation object,
      radecimage = TriScatteredInterp (alpha (sel), delta (sel), ... 
										real(skymap (sel)));

      % Create the regularly sampled RA/dec plane
      [ragrid, decgrid] = meshgrid (linspace (0,2*pi, uvsize), ...
									linspace (-pi/2,pi/2, uvsize));
          
      % generate samples from the interpolated model skyimage
      radecskymap = radecimage (ragrid, decgrid);
    else 
      % disp ('Not converting to RA/DEC images');
      radecskymap = lmskymap;
    end
