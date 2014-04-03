% Script to plot the track of A-team sources over a given timerange, for a zenith 
% pointing AARTFAAC.
% pep/17Mar14
% Arguments: 
%  starttime: UTC time of start of track (MJDsec)
%  endtime  : UTC time of end of track (MJDsec)
%        dt : Interval at which to generate Ateam positions (sec)
%        plt: bool controlling generation of plot of the Ateam track.
% Returns:
%  el : Elevation track of Ateam sources (rad).
%  az : Azimuth track of Ateam sources (rad).

function [el, az] = plotAteam (starttime, endtime, dt, plt)
	% Get the catalog positions of Ateam sources in equatorial coordinates
	load '~/WORK/AARTFAAC/Afaac_matlab_calib/srclist3CR.mat' % Load 3CR catalog.

	% AARTFAAC lat/long
    L = 6.869837540;       % longitude of CS002 in degrees
    B = 52.915122495;      % latitude of CS002 in degrees
	% Generate the list of times corresponding to the desired timerange and 
	% interval
	JD = [starttime:dt:endtime]./86400 + 2400000.5; % Convert MJDsec to JD.

	% Convert equatorial positions to azel.
	% First convert 
	a1 = 24110.54841;
	a2 = 8640184.812866;
	a3 = 0.093104;
	a4 = -6.2e-6;
	polcoeff = [a4, a3, a2, a1];
	
	% Time in Julian centuries
	TU = (floor(JD) + 0.5 - 2451545) / 36525;
	
	% Greenwich Star Time in seconds
	GST = (JD - floor(JD) - 0.5) * 86400 * 1.002737811906 + polyval(polcoeff, TU);
	
	% Local Siderderial Time in radians
	% 240: number of seconds per degree
	LST =  ((GST + L*240) / 240) * pi / 180;
	
	alpha0 = LST;              % Zenith projection, in radians
	delta0 = B * pi / 180;     % ibid.
	
	% Catalog Coordinates are B1950 radians
	% Offset of Ateam sources in 3CR catalog.
    ateam =  [324, 283, 88, 179]; 
	alpha = [srclist3CR(ateam).alpha];
	delta = [srclist3CR(ateam).delta];

	% Convert to cartesian
	[x, y, z] = sph2cart (alpha, delta, 1);  
	precMat = precessionMatrix(JulianDay(datenum(1950, 1, 1, 0, 0, 0)));	
	j2000cart = [x' y' z'] * precMat; % Transpose for product compatibility.
	% NOTE: 2*pi wraps in alpha do not matter.
	[alpha, delta, r] = cart2sph (j2000cart (:,1), j2000cart (:,2), ...
								  j2000cart (:,3));
	
	% (alpha, delta) to (l, m)
	el = zeros (length (ateam), length (alpha0));
	az = el;
	% For every time specified in the timerange:
	for ind = 1:length (alpha0)
		el(:, ind) = asin(sin(delta)*sin(delta0) + cos(delta0)*cos(delta) .* ...
						   cos(alpha0(ind) - alpha));
		az(:, ind) = acos((sin(delta)-sin(el(:,ind)) * sin(delta0)) ./ ...
						  (cos(el(:,ind)) * cos(delta0)));
		az(:, ind) = az(:, ind) .* (1 - 2 * (sin(alpha0(ind) - alpha) < 0));
	end;

	if (plt ~= 0)
%		plot (az(1,:)*180/pi, el(1,:)*180/pi, 'r'); hold on;
%		plot (az(2,:)*180/pi, el(2,:)*180/pi, 'b'); hold on;
%		plot (az(3,:)*180/pi, el(3,:)*180/pi, 'k'); hold on;
%		plot (az(4,:)*180/pi, el(4,:)*180/pi, 'g'); hold on;
		plot (az'*180/pi, el'*180/pi); 
		grid on;
		ylabel ('Elevation (deg)'); xlabel ('Azimuth(deg)');	
		title (sprintf ('Ateam tracks: %s to %s.',...
	    datestr (mjdsec2datenum (starttime)), datestr (mjdsec2datenum (endtime))));
		legend ('CasA', 'CygA', 'VirA', 'TauA');
		set(gca,'FontSize', 16,'fontWeight','bold')
		set(findall(gcf,'type','text'),'FontSize', 16, 'fontWeight','bold')
		

		% Plot azi and ele as a function of time
		dnum = mjdsec2datenum ([starttime:dt:endtime]);
		figure;
		subplot (211);
		plot (dnum, az*180/pi);
		grid on; axis tight;
		ylabel ('Azimuth (deg)'); xlabel ('Time');	
		datetick ('x', 15, 'keepticks');
		% title (sprintf ('Ateam tracks: %s to %s.',...
	    %    datestr (mjdsec2datenum (starttime)), datestr (mjdsec2datenum (endtime))));

		subplot (212);
		plot (dnum, el*180/pi);
		grid on; axis tight;
		ylabel ('Elevation (deg)'); xlabel ('Time');	
		datetick ('x', 15, 'keepticks');
		% title (sprintf ('Ateam tracks: %s to %s.',...
	    %    datestr (mjdsec2datenum (starttime)), datestr (mjdsec2datenum (endtime))));
		samexaxis ('join');
		p=mtit(sprintf ('Elevation and Azimuth tracks of Ateam sources between %s and %s.',...
	        datestr (mjdsec2datenum (starttime)), datestr (mjdsec2datenum (endtime))), ...
		   'xoff',-.07,'yoff',.015);
		legend ('CasA', 'CygA', 'VirA', 'TauA');
		set(gca,'FontSize', 16,'fontWeight','bold')
		set(findall(gcf,'type','text'),'FontSize', 16, 'fontWeight','bold')
	end;
