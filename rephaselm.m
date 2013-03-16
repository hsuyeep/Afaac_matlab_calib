% Script to rephase the passed ACM in the specified lm direction.
% pep/02Nov12
% Arguments:
%	acm  : ACM to rephase
%   l, m : Image plane coordinates to which to rephase
%   tobs : Time of observation in MJD sec. (for RA/Dec calculation)
%   freq : Frequency of observation in Hz.
% poslocal: Antenna positions in ITRF. Should be consistent with passed ACM.

function [re_acm, rot_uvw] = rephaselm (acm, l, m, tobs, freq, poslocal)
	tobs_jd = tobs/86400 + 2400000.5; % Convert MJD sec to JD

    % Generate uvw coords. in local coordinates (ie, wrt. CS002);
%	u = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
%	v = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
%	w = meshgrid (poslocal(:,3)) - meshgrid (poslocal(:,3)).';
%	uvw = [u(:) v(:) w(:)]; % 288*288 X 3 matrix.
%
%	% Convert passed l,m coords to theta (angle around u-axis) 
%	% and phi (angle around v-axis);
%	phi = acos (l); th = acos (m);

	% Convert passed l,m coords to RA/Dec.
	[r, d] = lmtoradec (l, m, tobs_jd, 6.869837540, 52.915122495); 
    % disp (['RA/dec of (l,m) = ' num2str([l,m]) ', is ' num2str([r d])]);

	% Convert RA/Dec to ITRF coords.
	rephase_pos = radectoITRF (r, d, false, tobs_jd); 

	% Generate array steering vector.
	A = exp (-(2*pi*1i*freq/299792458 * poslocal * rephase_pos.'));
	ph = A * A';

	% Steer the passed ACM in the chosen direction.
 	re_acm = ph .* acm;


%%%%%%%%%%%%%%  NOTE: One needs to rotate the uvw coordinates only when carrying
%%%%%%%%%%%%%% out facetted imaging, where the l,m coord. origin changes to
%%%%%%%%%%%%%% that of the center of the facet.
%	% Rotation matrix for rotation to reach the given l coord. of facet center,
%	% Via rotation about v-axis, by an angle corres. to phi.
%	rot_v = [[1        0         0]
%			 [0 cos(phi) -sin(phi)]
%			 [0 sin(phi)  cos(phi)]
%			];
%
%	% Rotation matrix for rotation about u-axis
%	% Rotation matrix for rotation to reach the given m coord. of facet center,
%	% Via rotation about u-axis, by an angle corres. to th.
%	rot_u = [[cos(th)  0  sin(th)]
%			 [0        1        0]
%			 [-sin(th) 0  cos(th)]
%			];
%
%	% Rotate the uvw coordinates of each baseline.
%	% NOTE: Can this be done by using the HA/Dec. relation to u,v,w coords., wit
%	% antenna positions in ITRF?
%	rot_uvw = rot_v * rot_u * uvw;
