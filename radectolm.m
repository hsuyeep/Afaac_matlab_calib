function [l, m] = radectolm(alpha, delta, JD, L, B)

% [l, m] = radectolm(alpha, delta, JD, L, B)
%
% Converts (ra, dec) coordinates to (l, m) coordinates assuming the
% specified time and the geographical location of Dwingeloo.
%
% Arguments
% alpha : vector or matrix with RA coordinates in radians
% delta : vector or matrix with the same size as alpha (not checked by the
%         program) with corresponding DEC coordinates in radians
% JD    : Julian day as DATENUM for the time of observation
% L     : geographical longitude in degrees
% B     : geographical latitude in degrees
%
% Return values
% l : vector or matrix of the same size as alpha with l-coordinates
% m : vector or matrix of the same size as alpha with corresponding
%     m-coordinates
%
% SJW, 2003
%
% SJW, 14 June 2011: Updated GST calculation to make it identical to the
% calculation in radectoITRF (only cosmetic, difference less than 0.01s).

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

% (alpha, delta) to (l, m)
el = asin(sin(delta) * sin(delta0) + cos(delta0) * cos(delta) .* cos(alpha0 - alpha));
az = acos((sin(delta) - sin(el) * sin(delta0)) ./ (cos(el) * cos(delta0)));
az = az .* (1 - 2 * (sin(alpha0 - alpha) < 0));
l = -real(cos(el) .* sin(az));
m = real(cos(el) .* cos(az));
l(el < 0) = NaN;
m(el < 0) = NaN;
