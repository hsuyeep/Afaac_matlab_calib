function [ra, dec] = SunRaDec(JD)

% [ra, dec] = SunRaDec(JD)
%
% Calculate the (ra, dec) position of the Sun in J2000 coordinates. The
% implementation is based on the Astronomical Almanac 2012.
%
% Argument
% JD  : time of observation in Julian days
%
% Return values
% ra  : right ascension of the Sun in radian
% dec : declination of the Sun in radian
%
% SJW, 22 June 2011
% The code has been verified by comparing its results to the positions
% listed in the Astronomical Almanac.

% number of days since J2000
n = JD - 2451545.0;
% mean longitude of the Sun, corrected for aberration, in degrees
L = mod(280.460 + 0.9856474 * n, 360);
% mean anomaly of the Sun, in degrees
g = mod(357.528 + 0.9856003 * n, 360);
% ecliptic longitude, in degrees
lambda = L + 1.915 * sind(g) + 0.020 * sind(2 * g);
% obliquity of ecliptic, in degrees
epsilon = 23.439 - 0.0000004 * n;
% right ascension, J2000, in degrees
ra = atan2(cosd(epsilon) * sind(lambda), cosd(lambda));
% declination, J2000, in degrees
dec = asin(sind(epsilon) * sind(lambda));
