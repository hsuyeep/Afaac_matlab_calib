% Script to convert time in UTC to MJD seconds
% Taken from pg. 7 of 'Practical Astronomy with your calculator'
% By Peter Duffet-Smith
% NOTE: The output does not match example taken from a measurement set!
% 2011-09-21-12:39:08 = 4823325543.237480
% pep/07Nov12

function [mjdsec, jd] = utc2mjdsec (y, m, d, hr, mn, sec)
	if m == 1 | m == 2
		y = y - 1;
		m = m + 12;
	end;
	a = double(int32(y/100));
	b = 2 - a + double (int32 (a/4));

	if y < 0
		c = double (int32 (365.25*y - 0.75));
	else
		c = double (int32 (365.25*y));
	end;
	e = double (int32(30.6001 * (m+1)));
	jd = (b + c + e + d + 1720994.5);
	mjdsec = (jd  - 2400000.5) * 86400 + hr*3600 + mn*60+sec;
