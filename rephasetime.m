% Program to rephase a given ACM towards a specified time offset.
% NOTE: Since earth rotation is only along u-axis in the visibility plane,
% the rephasing needs to be done only in 1 dimension.
%   based on the fourier transform shift theorem:
%   	f(x-dx,y-dy) <==> exp(-2pi i(u*dx+v*dy)) F(u,v)
%
%
% Arguments:
%	acm : The ACM to rephase
%	tobs: Time of observation of the acm, in MJD seconds.
%	trephase: Time, in MJD seconds, to which to rephase.
%	freq: Frequency of observation of the acm, in Hz.
%	posITRF: ITRF Antenna positions.
%
% Returns :
%	rephase_acm: The acm after rephasing to the desired time.
%
% NOTE: Typically, each ACM is integrated over one second. For the LBA, at the 
% center frequency of 60 MHz (5m wavelength), there exists a 10deg. phase diff.
% for a source ~ TODO
% 
function [re_acm] = rephasetime (acm, tobs, trephase, freq, posITRF)
	c = 299792458; % m/s
	tobs_jd = tobs/86400 + 2400000.5;
	treph_jd = trephase/86400 + 2400000.5;

	% Determine RA/DEC of zenith (this obs)
	[r, d] = lmtoradec (0, 0, tobs_jd, 6.869837540, 52.915122495); 

	% Determine RA/DEC of point to rephase to..
	r = r + (tobs - trephase) * (2*pi/86400); % rad

	% Convert to ITRF 
	rephase_pos = radectoITRF (r, d, false, tobs_jd); 

	% Generate array steering vector.
	A = exp (-(2*pi*1i*freq/299792458 * posITRF * rephase_pos.'));
	ph = A * A';

	% Steer the passed ACM in the chosen direction.
 	re_acm = ph .* acm;

%	% Determine the zenith angle (in radians) corresponding to trephase.
%	za = ((trephase-tobs)*2*pi)/86400; 
%	ph = 2*pi*freq/c * sin(za);
	% x_shift = exp(-i*2 * pi * delta(1) * [0:floor(N/2)-1 floor(-N/2):-1]' / N);
	% y_shift = exp(-i*2 * pi * delta(2) * [0:floor(M/2)-1 floor(-M/2):-1] / M);

	% Form the phase matrix towards the direction specified by trephase
	% for every baseline. The sign of the phase is correctly taken from the
	% position difference (poslocal).
%	uloc = ph * meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
%    % vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
%	rephase_acm = uloc .* acm;
