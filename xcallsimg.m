% Script to extract out the flux of chosen sources using LS imaging on 
% calibrated snapshot visibilities.
%	The main reason for having a separate function for LS imaging, and not 
% accomodating LS imaging in extractcal.m is its operation on calibrated vis.
% pep/17Jul13
%  Arguments:
%	acc : The calibrated array correlation matrix
%	tobs: Time of observation in MJD secs
%   freq: Freq. of observation in Hz.
% rodata: Structure containing flagged posITRF positions (consistent with acc).
%catalog: Catalog to use for positional information of sources.
%callist: List of monitor sources.
% epoch : True implies catalog positions are B1950
% debug : Flag turning on the display of debug information.
% Returns :
%	upcals: Selection vector for list of calibrator sources above the horizon.
%   fit : Structure containing estimated flux information for each source.

function [upcals, fit] = xcallsimg (acc, tobs, freq, rodata, catalog,...
									 callist, epoch, debug)
    tobs_jd = tobs/86400 + 2400000.5; % Convert from MJD secs. to JD day units
	% Extract out calibrators within FoV and above horizon
	racal  = [catalog(callist).alpha];
	deccal = [catalog(callist).delta];
	fluxcal= [catalog(callist).flux];

	srcposITRF = radectoITRF (racal, deccal, epoch, tobs_jd);
	upcals = srcposITRF * rodata.normal > 0;
	% Generate the array steering vector for all selected sources.
	A = exp(-(2 * pi * 1i * freq / 299792458) * ... 
				(rodata.posITRF_fl * srcposITRF (upcals, :).'));

	[catl, catm] = radectolm (racal,deccal,tobs_jd,rodata.lon,rodata.lat,epoch);
	fprintf (1, '%.2fs %.2fHz: Found %d calibrators in FoV.\n', ...
			 tobs, freq, sum(upcals));
	cat = [catalog(callist(upcals))];
	lup = catl(upcals); mup = catm(upcals);

	if debug > 0
		for ind = 1:sum(upcals)
			fprintf (1, 'Name: %s, flux: %.2f Jy, l: %8.5f, m: %8.5f\n', ...
					cat(ind).name, cat(ind).flux, lup(ind), mup(ind));
		end;
	end;

	% Estimate the flux of each source via least squares imaging
	numnans= sum (sum (isnan(acc))); % List the number of usable baselines.
	acc (isnan(acc)) = 0; % Convert NaNs to 0s before using acc;
	fit.peakfl = numnans * real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * acc (:));
