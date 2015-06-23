% Script to attempt spherical harmonic imaging from AARTFAAC.
% Arguments:
%	simparm : Array simulation parameters, used to call simvis(). See that function for parameters;
%	deb     : Bool to turn on debugging plots.
% Based on [1]: Carozzi, 'Imaging on a sphere with interferometers: the spherical wave harmonic transform'
%				arxiv:1504.04485v2, Apr. 2015
% Returns:
%	swht  : The harmonic components of the visibilities.
%	blm   : The brightness distribution's spheric harmonic components.
%	img   : The reconstructed image from the spherical harmonic components.
% pep/21May15

% NOTE: Sample to plot polar functions in matlab
% azi = linspace (0, 360, 361)*pi/180; % Uniform sampling of azi and ele. axis
% ele = linspace (0, 180, 181)*pi/180;
% [A, E] = meshgrid (azi, ele);
% x_deb = sin(A).*cos(E); % Convert to cartesian coordinates.
% y_deb = cos(A).*cos(E);
% z_deb = sin(E);
% plot3(x_deb(:), y_deb(:), z_deb(:), '.')
% surf (x_deb, y_deb, z_deb, ones (size(x_deb))); % Colormap intentionally set to 1
% shading flat

function [swht, blm, img] = swhtimg(simparm, deb)

	if (isempty (simparm))
		fname = '~/WORK/AARTFAAC/Reobs/20Nov13/r02/SB002_LBA_OUTER_8b2sbr02_1ch.bin';
		load ('poslocal_outer.mat', 'poslocal');
	end;
	doplot = 0;
	
	Nelem = 288;
	fid = fopen (fname, 'rb');
	[acc, t_obs, freq] = readms2float (fid, 100, -1, Nelem);

	 % Wavenumber for this observation.
	lambda = 299792458./freq;
	k0 = (2*pi)/lambda; % 1m in units of phase (rad) of a wave of wavelength lambda.
	
	lmax = 4; % int32(Nelem/50); % Highest band of harmonic functions to be implemented.
	% Store for coefficients of the SWH ananlysis. +1 due to lmax going from 0-lmax.
	% NOTE: Coefficients are stored in the lower triangle of the swht matrix. 
	% Each row = different band (l-value). First column is m=-l, till m=l
	swht = zeros (lmax+1, 2*lmax+1); 
	
	% Choose set of antennas to work with. Currently 1 station for ease of imaging.
	stationsel = [1:48];

	% Baseline coordinates.
	uloc = meshgrid (poslocal(stationsel,1)) - meshgrid (poslocal(stationsel,1)).';
	vloc = meshgrid (poslocal(stationsel,2)) - meshgrid (poslocal(stationsel,2)).';
	wloc = meshgrid (poslocal(stationsel,3)) - meshgrid (poslocal(stationsel,3)).';
	
	% Convert cartesian coordinates of each baseline to spherical coordinates.
	[blth, blphi, blR] =  cart2sph (uloc(:), vloc(:), wloc(:)); % ITRF m, local to CS002.
	blR = blR*freq/299792458; % Radial component in lambda units.
	blR = blR * k0; % Baseline radial extent, expressed as phase (in radians).
	
	% ACM with selected stations.
	acc = acc (stationsel, stationsel);
	if (doplot == 1)
		imagesc (abs(acc)-diag(diag(acc))); colorbar;
		title (sprintf ('ACM for time %.0f\n', tobs));
	end;

	% Convert to a vector
	acc = acc (:);
	
	% Initialize stuff for debug plots
	if (deb > 0)
		azi = linspace(0, 360, 361) * pi / 180;
		ele = linspace(0, 180, 181) * pi / 180;
		[A, E] = ndgrid(azi, ele);
        % Observed spherical coordinates of the visibilities converted to
        % cartesian
        x_deb = cos(blphi).*sin(blth);
        y_deb = sin(blphi).*sin(blth);
        z_deb = cos(blth);        
	end;
	% Compute the spherical harmonic transform coefficients over the chosen quantal number range.
	for l = 0:lmax
		fprintf (2, '.');
		for m = -l:l
            fprintf (1, '  (l=%d, m=%d)\n', l, m);
			mabs = abs(m); % Since spherical harmonics for the m-quantal numbers are symmetric.

			% Analyze each visibility into spherical harmonics. 
			% Each SH coefficient is a weighted sum of all visibilities.
			% c.f. eq. 16 of [1]
			for ind = 1:length(acc)             
				% Generate the spherical bessel function for these quantal numbers.
				if (blR(ind) == 0) % For the zero spacing case
					if (l == 0) 
						bes = 1;
					else
						bes = 0;
					end;
				else
					% NOTE: This is the spherical bessel function of the first kind, 
					% sampled at the radial extent of the baseline.
					% bes = sqrt (pi./(2*blR(ind)*k0)).*besselj(l+0.5, blR(ind)*k0);

					% Since radial component of the baseline is already in cycle units.
					bes = sqrt (pi./(2*blR(ind))).*besselj(l+0.5, blR(ind)); 
				end;
	
				% Generate the spherical harmonic coefficients for these quantal numbers
				Plm = legendre(l,cos(blth(ind)));
				if l~=0
					Plm = squeeze(Plm(mabs+1,:,:));
				end
				a1 = ((2*l+1)/(4*pi));
				a2 = factorial(l-mabs)/factorial(l+mabs);
				C = sqrt(a1*a2);
				Ylm = C*Plm.*exp(i*mabs*blphi(ind));
				if (m < 0)
					Ylm = conj(Ylm)*(-1)^mabs;	
				end;
				
				if (deb > 0)
					Ylmdeb(ind) = Ylm;
				end;
				swht (l+1,m+l+1) = swht (l+1,m+l+1) + acc(ind)*bes*Ylm;
			end;
            if (deb > 0)
                subplot (121);
                plot (blth, real(Ylmdeb), '.');
                xlabel ('baseline theta'); ylabel ('Real Ylm');
                title (sprintf ('Y(l,m): l=%d, m=%d',l,m));
               
                subplot (122);
                plot (blth, imag(Ylmdeb), '.');
                xlabel ('baseline theta'); ylabel ('Imag Ylm');
                title (sprintf ('Y(l,m): l=%d, m=%d',l,m));
             
                % Convert the sampled spherical harmonic into cartesian coords
                Data = abs(imag(Ylmdeb));
                minData = min(Data(:));
                maxData = max(Data(:));
                Distord = (Data - minData)/(maxData-minData);
                % surf (x_deb, y_deb, z_deb, Data);
                % shading flat;
                % title (sprintf ('Y(%2d, %2d)',l, m));
            end;
		end;
	end;
	
	fprintf (1, '\n<-- Computed shwt coe.\n');
	if (deb > 0)
		imagesc ([0:lmax], [-lmax:lmax], abs(swht)); colorbar;
		xlabel ('m-quantal'); ylabel ('l-quantal');
		title ('SWH coefficients for snapshot visibility.\n');
	end;

	% Compute the spherical harmonic coefficients of the brightness distribution on the sky.
	% c.f. eq. 11 of [1]
	for ind = 1:size (swht,1)
		blm(ind, :) = swht(ind, :)/(4*pi)*(-i)^(-ind);
	end;
	
	% Reproject to local coordinates for generating an image in local coordinates.
	res = 0.05; % Resolution of image in direction cosine units.
	lloc = meshgrid ([-1:res:1]);
	mloc = lloc';

	% Convert the 2-D cartesian coordinates of the local image to 2-D polar
	% to represent the image as a sum of spherical harmonics.
	[phisph, rsph] = cart2pol(lloc, mloc);
	thsph = asin (rsph);
	% thsph = acos (rsph);
	img = zeros (size (lloc));
	
	for xpix = 1:size (img,1)
		fprintf (2, '+');
		for ypix = 1:size (img,2)
			if (rsph(xpix, ypix) <= 1) % Only image those pixels within the circle in local coords.
				for l = 0:lmax
					for m = -l:l
						mabs = abs(m);
						% Generate the spherical harmonic coefficients for these quantal numbers
						Plm = legendre(l,cos(thsph(xpix,ypix)));
						if l~=0
							Plm = squeeze(Plm(mabs+1,:,:));
						end
						a1 = ((2*l+1)/(4*pi));
						a2 = factorial(l-mabs)/factorial(l+mabs);
						C = sqrt(a1*a2);
						Ylm = C*Plm.*exp(i*mabs*phisph(xpix,ypix));
						if (m < 0)
							Ylm = conj(Ylm)*(-1)^mabs;	
						end;
				
						% c.f. eq. 6 of [1]
						img(xpix, ypix) = img(xpix, ypix) + blm(l+1,m+l+1)*Ylm;	
					end;
				end;
			end;
		end;
	end;
	fprintf (1, '\n<-- Done imaging with SPH coefficients.\n');


