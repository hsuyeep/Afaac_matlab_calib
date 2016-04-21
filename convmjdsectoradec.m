% Script to convert the position in the sky corresponding to the zenith location
% for the Superterp.
% Based on the approximate method in http://stargazing.net/kepler/altaz.html
% LST = 100.46 + 0.985647 * d + long + 15*UT
%      d    is the days from J2000, including the fraction of
%           a day
%      UT   is the universal time in decimal hours
%      long is your longitude in decimal degrees, East positive.
%      
% Add or subtract multiples of 360 to bring LST in range 0 to 360
% degrees.
% Arguments:
%	mjdsec : Time of observation in MJD seconds.
% Returns  :
%  ra,dec  : RA/Dec of zenith in radians.
% pep/07Nov12

function [ra, dec] = convmjdsectoradec (mjdsec)
	% location of array
	lon     = 6.869837540;       % longitude of CS002 in degrees, East = +ve 
    lat     = 52.915122495;      % latitude of CS002 in degrees, North = +ve

	% Days from J2000 corresponding to obs time in MJD sec
	% J2000 = JD 2451545.0
	jd2j2000 = int32((mjdsec/86400.) - (2451545.0-2400000.5)); 
	jd2j2000 = double (jd2j2000); % Else all calculations are also int. only!
	% Fraction of day remaining, in UTC
	fracday = ((mjdsec/86400.) - (2451545.0-2400000.5)) - jd2j2000; 

	lst = 100.46 + 0.985647 * jd2j2000 + lon + 15*fracday*24; % In degrees
	% Eliminate multiples of 360deg. to bring within 2pi
	lstturns = double(int32(lst/360)); 
	lst = lst - lstturns * 360;
	ra  = lst * pi / 180 + pi; % For HA = 0 (zenith), RA = LST.
	dec = lat * pi / 180; % Dec of zenith = lat of array.
