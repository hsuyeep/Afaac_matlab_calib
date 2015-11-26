% Script to convert a given mjdsec to unixtime
% pep/25Nov15

function utime = mjdsec2unixtime (mjdsec)
	utime = (mjdsec2datenum (mjdsec) - datenum ([1970 1 1]))*86400;
