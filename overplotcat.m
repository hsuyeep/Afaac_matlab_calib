% Script to overplot 3CR catalog sources over the current image. 
% Returns a list of sources with offsets in catalog in sel, and positions in 
% sel_l, sel_m.
% pep/24Oct12
% NOTE: the 3CR catalog (srclist3CR.mat) contains sources positions in B1950 
% coordinates, and in radians, while the VLSS catalog contains positions in 
% J2000 coordinates. Hence an appropriate conversion needs to be done.

% Arguments:
%	tobs : Time of observation, in MJD secs.
% catalog: The 3CR catalog of sources.
% abovejy: The flux cutoff, in Jy@158MHz for the number of sources.
% hdl    : Handle of the plot in which to overplot.
% epoch  : bool: true for catalog positions in B1950, false for J2000 positions.

% Returns :
%	sel_[l,m]: The l,m coordinates for selected sources.
% 	sel      : The catalog offsets of the selected sources (to get other info. 
%			   from the catalog, e.g. names and flux.

function [sel, sel_l, sel_m] = overplotcat (tobs, catalog, abovejy, hdl, epoch)
	% Overplotting 3CR sources above 'abovejy' Jy, in l,m coordinates
	% hold on;
	tobs_jd = tobs/86400 + 2400000.5; % Convert from MJD secs. to MJD day units
	cs002_lon     = 6.869837540;      % longitude of CS002 in degrees
    cs002_lat     = 52.915122495;     % latitude of CS002 in degrees
	cs002_normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	
	% Choose sources above flux cut off limit.
	sel = ([catalog(:).flux] > abovejy);
	
	% Choose sources above current horizon. 
	srcpos0 = radectoITRF ([catalog(sel).alpha], [catalog(sel).delta], ... 
							epoch, tobs_jd);
    up = srcpos0 * cs002_normal > 0;
	sel (sel == 1) = up;

	[sel_l, sel_m] = radectolm ([catalog(sel).alpha], ...
		[catalog(sel).delta], tobs_jd, cs002_lon, cs002_lat, epoch);

	% plot (srcm_25, srcl_25, 'o'); 
	namedsrc = [catalog(sel)];
	
	% Take care of l.m zenith projection by ignoring sources at the edges.
	for i = 1:length (namedsrc)
	  if namedsrc(i).flux > abovejy && abs(sel_m(i))<0.99 && abs(sel_l(i))<0.99
	    % text ('Color' , [1 1 1], 'Position', [sel_m(i) sel_l(i)], ... 
	    text ('Color' , [1 1 1], 'Position', [sel_l(i) sel_m(i)], ... 
			  'String', namedsrc(i).name);
	  end
	end	
