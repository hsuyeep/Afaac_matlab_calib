% Script to convert a UTC time in datenum to MJDsec, conforming to LOFAR usage.
% pep/05Nov15

% Arguments:
%   tobs  : UTC time as a datenum

% Returns:
%   mjdsec : Specified time as MJDsec

function mjdsec = datenum2mjdsec (tobs)
	% Convert MJD reference date to datenum
	mjddateref = datenum (1858,11,17,00,00,00); % Start of MJD
    mjdsec = (tobs - mjddateref)*86400.;
