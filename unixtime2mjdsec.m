% Script to convert unixtime to MJDsec
% pep/06Aug15
% Arguments:
%	unixtime: Time to be converted, specified in seconds past 01Jan1970.
% Returns:
%	mjdsec: Specified time in mjdsec seconds.

function mjdsec = unixtime2mjdsec (unixtime)
	mjddateref = datenum (1858,11,17,00,00,00); % Start of MJD
	mjdsec = (unixtime2datenum(unixtime) - mjddateref)*86400;
