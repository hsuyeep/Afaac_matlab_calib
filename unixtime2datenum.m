% Script to convert unix time (seconds past 01Jan1970) to datenum
% pep/06Aug15
% Arguments:
%  unixtime: Time to convert in seconds past 01 Jan 1970
% Returns:
%  dnum : Matlab time corresponding to the passed unixtime.

function dnum = unixtime2datenum (unixtime)
	dnum = unixtime/86400. + datenum (1970, 1, 1);
