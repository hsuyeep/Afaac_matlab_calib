% Script to extract the RA of the zenith from LOFAR CS002, corresponding 
% to a given time in UTC.
% pep/24Jul15
% Arguments:
%   utc : String containing the UTC time in the form of an array. [YYYY, MM, DD, HH, MM, SS]

function [LST, JD] = utc2ra (utc)
	CS002_lon     = 6.869837540;       % longitude of CS002 in degrees
    CS002_lat     = 52.915122495;      % latitude of CS002 in degrees

	utcdatenum = datenum (utc);
	JD = utcdatenum - datenum(1998, 2, 13, 12,0,0) + 2450858;

	a1 = 24110.54841;
	a2 = 8640184.812866;
	a3 = 0.093104;
	a4 = -6.2e-6;
	
	% Code taken from radectolm.m for calculation of LST.
	polcoeff = [a4, a3, a2, a1];
	
	% Time in Julian centuries
	TU = (floor(JD) + 0.5 - 2451545) / 36525;
	
	% Greenwich Star Time in seconds
	GST = (JD - floor(JD) - 0.5) * 86400 * 1.002737811906 + polyval(polcoeff, TU);
	
	% Local Siderderial Time in radians
	% 240: number of seconds per degree
	LST =  ((GST + CS002_lon*240) / 240) * pi / 180;
	LST = mod(LST, 2 * pi);
	fprintf (1, 'JD (days): %.4f, LST (rad): %.4f\n', JD, LST);

