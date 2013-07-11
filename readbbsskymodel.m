% Script to read in a BBS Solar sky model (a text file) of the Sun. The 
% position of the Sun during generation of the Sky model is fixed (hardcoded).
% The code translates the Solar model to the position of the Sun specified by 
% tobs, and generates ITRF positions of the components, along with a 
% normalized flux. In case the model components are gaussians, the gaussian 
% parameters are returned.

% The format string parsed is:
% format = Name, Type, Ra, Dec, I, Q, U, V, MajorAxis, MinorAxis, Orientation, ReferenceFrequency='1.50397e+08', SpectralIndex='[]'
% The components handled are of type POINT and GAUSSIAN.
% pep/29May13

% Arguments:
%  fname : Filename of text file containing src model components
%   tobs ; Time of observation in MJD.
%    hba : Bool indicating the presence of an HBA sky model in fname. else LBA
%  debug : Bool for turning on debugging, displays components as a plot.

% Returns:
%  srcitrf: An array of structures containing the parameters of the model 
%           components, with the position in ITRF.

function [srcitrf] = readbbsskymodel (fname, tobs, hba, debug);

	% Positional information.
	cs002_lon     = 6.869837540;      % longitude of CS002 in degrees
    cs002_lat     = 52.915122495;     % latitude of CS002 in degrees
    cs002_normal  = [0.598753, 0.072099, 0.797682].'; % Unit normal in ITRF to CS002

	% NOTE: datenum returns seconds past 01Jan0000.
	% Taking center of time window
	if (hba == 0)
		% 23Jun2012, 11:04:52.1-11:05:20.4 UTC, LBA
		% per email of 04Jun13 by Vocks to p.prasad@uva.nl.
		% tobs_mod = JulianDay (datenum ([2012, 6, 23, 11, 5, 6.25]));
		tobs_mod = JulianDay (datenum ([2012, 6, 23, 11, 4, 52]));
		tit_str = 'LBA 9 component Solar model, 23/06/2012, 11:05:6.25UT';
	else
		% 23Jun2012, 11:19:42.6 - 11:20:11.0 UTC in HBA
		tobs_mod = JulianDay (datenum ([2012, 6, 23, 11, 19, 56.8]));
		tit_str = 'HBA 34 component Solar model, 23/06/2012, 11:20:11UT';
	end;

	if (isempty (fname))
		error ('Empty file!');
	end;

	fid = fopen (fname, 'rt');
	
	% Read in the file header
	head = textscan (fid, 'format = %s %s %s %s %s %s %s %s %s %s %s %s=%f SpectralIndex=''[]''', 1);
	fprintf (1, 'Reference frequency: %s\n', char(head{12}));

	% Read in blank line after the header.
	fgets (fid);

	% Read in all components
	comp = textscan (fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s');
	cmps = length (comp{1});

	fprintf (1, 'Number of components found: %d. Fields: %d\n', ... 
			 cmps, length (comp));

	% Find current mean position of the Sun, as well as during the model sky
	% generation. NOTE: Returned RA/DEC are in  Units of radians.
    tobs_jd = tobs/86400 + 2400000.5; % Convert from MJD secs. to JD days
	[raSun, decSun] = SunRaDec (tobs_jd);
	[raModSun, decModSun] = SunRaDec (tobs_mod);
	raModSun  = raModSun * 180/pi*3600; % Convert to arcsec
	decModSun = decModSun * 180/pi*3600; % Convert to arcsec
	raSun     = raSun * 180/pi*3600; % Convert to arcsec
	decSun    = decSun * 180/pi*3600; % Convert to arcsec

	if (debug == 1)
		fprintf (1, 'Model Sun: %f, %f, Current Sun %f, %f\n', ...
				 raModSun/3600, decModSun/3600, ...
				 raSun/3600, decSun/3600);
	end;
	
	
	% Convert components to decimal RA/DEC.
	% Assumed format: HH:MM:sec.sec for RA, DEG.min.sec.sec for DEC.
	ra = zeros (cmps, 1); dec = ra;
	for cmp = 1:cmps
		[str, remra] = strtok (comp{3}(cmp), ':');	% Get hours
		ra(cmp) = str2num (char(str));
		[str, remdec] = strtok (comp{4}(cmp), '.');
		dec(cmp) = str2num (char(str));
		for ind = 1:2 
			[str, remra] = strtok (remra, ':');	
			ra(cmp) = ra(cmp) + str2num (char(str))/(60^ind);
			[str, remdec] = strtok (remdec, '.');
			dec(cmp) = dec(cmp) + str2num (char(str))/(60^ind);
		end;
	end;
	fprintf (2, 'NOTE: Fractional secs in dec not considered yet!\n');

	if (debug == 1)
		for ind = 1:length (ra)
			fprintf (1, 'Comp %d: %f, %f\n', ind, ra(ind), dec(ind));
		end;
	end;

	% Convert RA from hours to arcsec
	ra = ra*(2*pi)/24*180/pi*3600;		 
	% Convert dec from degrees to arcsec
	dec = dec * 3600; 

	% Move model components to the position of the Sun in our observations.
	ra_curr = ra - raModSun + raSun;
	dec_curr = dec - decModSun + decSun;

	% Convert the positions of all components to ITRF
	% pos = radectoITRF (ra, dec, false (length(dec), 1), tobs_jd);

	% Convert the positions of all components from RA/DEC to l,m for a 
	% zenith pointing from CS002 at the time specified in tobs.
%	pos = zeros (length (ra), 2);
%	[pos(:,1), pos(:,2)] = radectolm (ra, dec, tobs_mod, cs002_lon, cs002_lat, false);
%	[sunmodl, sunmodm] = radectolm (raModSun, decModSun, tobs_mod, cs002_lon, cs002_lat, false);
%
%	[sunl, sunm] = radectolm (raSun, decSun, tobs_jd, cs002_lon, cs002_lat, false);

	% Obtain I flux, normalize to 1, scale by arbit factor for plot.
	iflux = str2num(char(comp{5}));
	normflux = (iflux ./ max (iflux))*1000;

	if debug == 1
		extent = 300; 
		% Extract out the gaussian parameters of each component.
		% List of all gaussian components
		tf = strcmp ('GAUSSIAN,', comp{2});

		% Major FWHM in arcsec
		maj =  (str2num(char(comp{9}(tf))));
		% Minor FWHM in arcsec
		minor = (str2num(char(comp{10}(tf))));

		% Orientation of Major axis in deg, 'North over East' possibly 
		% implies North is 0, and inc. is in the East direction(clockwise).
		pa = pi/180*str2num(char(comp{11}(tf))); % Convert to rad.

		% Determine the l,m widths of major and minor axis.
		% TODO

		% Generate gaussians with specified parameters in the lm plane.
		a = cos(pa).^2./(2*maj.^2) + sin(pa).^2./(2*minor.^2);
		b = -sin(2*pa)./(4*maj.^2) + sin(2*pa)./(4*minor.^2);
		c = sin(pa).^2./(2*maj.^2) + cos(pa).^2./(2*minor.^2);
		
		% Form a grid in RA/DEC space. The 2D gaussians get sampled on this.
		% rares = 0.001; decres = 0.001; % In radians.
		rares = 10; decres = 10; % In arcsec
		% NOTE: Going 300 arcsec beyond available components.
		x = min(ra)-300:rares:max(ra)+300;  
		y = min(dec)-300:decres:max(dec)+300; 

		% Form the grid on which sources will lie.
		[xg, yg] = meshgrid (x, y);

		% Stores the model with all components.
		model = zeros (size (xg));

		% Stores individual components, imaged over same extent as model
		gaucomps = zeros (size (xg,1), size(xg,2), cmps);

		for ind = 1:cmps
			gaucomps (:,:,ind) = normflux (ind) * exp (-(a(ind)*(xg-ra(ind)).^2 + 2*b(ind)*(xg-ra(ind)).*(yg-dec(ind)) + c(ind)*(yg-dec(ind)).^2));
			model = model + gaucomps (:,:,ind);
		end;

		% Generate actual solar disk of 0.5deg
		sunmask = ones (size (xg));
		sunmask(((xg-raModSun).^2+(yg-decModSun).^2)<...
			770^2 & ...
		((xg-raModSun).^2+(yg-decModSun).^2)>750^2)=0;

		% Plot the composite model, add in the model Sun position
		model = model .* sunmask;
		sunmask (sunmask == 0) = 200;
		model = model + sunmask;
		figure;
		imagesc (x/3600, y/3600, model + sunmask); colorbar;
		hold on;
		plot (raModSun/3600, decModSun/3600, 'r*');
		title (strcat (tit_str,' :Composite'));
		xlabel ('RA (deg)');
		ylabel ('DEC(deg)');

		% Plot the individual components with chosen extent
		model = zeros (size (xg));
		figure;
		for ind = 1:cmps
			mask = zeros (size (xg));
			mask (((xg-ra(ind)).^2 + (yg-dec(ind)).^2) < extent^2) = 1;
			model = model + gaucomps (:,:,ind).*mask;
		end;
		imagesc (x/3600, y/3600, model + sunmask); colorbar;
		hold on;
		plot (raModSun/3600, decModSun/3600, 'r*');
		title (tit_str);
		xlabel ('RA (deg)');
		ylabel ('DEC(deg)');
	end;
