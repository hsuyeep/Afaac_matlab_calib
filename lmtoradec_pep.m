% Script to convert provided l,m coordinates to RA/dec around zenith
% pep/12Mar12

function [ra, dec] = lmtoradec (l, m, JD, L, B)
% l,m = vector of l and m positions, both of equal size
% JD  = JD of exact time of observation, in days
% L   = longitude of observing position, in decimal degrees
% B   = latitude of observing position, in decimal degrees

%% Locate the RA/dec of the zenith point
deg2rad = (pi/180); % Multiply by this to convert a value in deg., into radians
% LST calculation (from radectolm.m)
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

alpha0 = LST - 2*pi*(floor (LST/(2*pi))); % Zenith projection, in radians between [0,2*pi]


% alpha0 = (100.46 + 0.985647 * JD + L)*deg2rad; % RA of zenith (radians) = LST of observation
delta0 = B*deg2rad;                % declination of zenith (radians) = latitude of observing location
lrad = (pi/2) * l;        % Convert l and m into radians, so that they range from -pi/2 to pi/2
mrad = (pi/2) * m;
sel  = ((1-lrad.*lrad - mrad.*mrad) > 0); % To avoid generating complex numbers while taking the sqrt.

%% Calculate the RA and dec for every l,m point as offsets to zenith point
ra = alpha0 + atan2(lrad, sqrt (1-lrad.*lrad - mrad.*mrad) .* cos (delta0) - mrad.*sin (delta0));
dec = asin (mrad .*cos(delta0) + sqrt (1 - lrad.*lrad - mrad.*mrad) .* sin(delta0));
