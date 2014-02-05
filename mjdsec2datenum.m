% Script to convert time in MJDseconds to datenum format.
% pep/24Jan14
% Arguments:
%	tobs: Time to be converted, specified in MJD seconds.
% Returns:
%	mjddatenum: Specified time in datenum format.

function mjddatenum = mjdsec2datenum (tobs)
	% Convert MJD reference date to datenum
	mjddateref = datenum (1858,11,17,00,00,00); % Start of MJD
	mjddatenum = mjddateref + tobs/86400.;
