% Script to extract out the flux of chosen sources using LS imaging on 
% calibrated snapshot visibilities, after positions are estimated using WSF. 
% Functional equivalent of extractcal.m.
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

	% Extract out monitor sources within FoV and above horizon
	racal  = [catalog(callist).alpha];
	deccal = [catalog(callist).delta];
	fluxcal= [catalog(callist).flux];

	srcposITRF = radectoITRF (racal, deccal, epoch, tobs_jd);
	upcals = srcposITRF * rodata.normal > 0;

	% Convert positions from ITRF to azi/el
	thsrc0 = asin (srcposITRF(upcals, 3));
	phisrc0 = atan2 (srcposITRF(upcals, 2), srcposITRF(upcals, 1));
	
	% Estimate positions of all sources using WSF, using unit gains, 0 tsys.
	[thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, fval, out] = ... 
    	wsf_srcpos (acc, ones (1,length(acc)), freq, zeros (size (acc)), ...
			sum(upcals), rodata, phisrc0, thsrc0, upcals, [], 1);

   % Convert WSF positions from azi/el to ITRF
    srcposhat = [cos(phisrc_wsf(upcals)) .* cos(thsrc_wsf(upcals)),...
                 sin(phisrc_wsf(upcals)) .* cos(thsrc_wsf(upcals)),...
                 sin(thsrc_wsf(upcals))];

	% Generate the array steering vector for all selected sources.
	A = exp(-(2 * pi * 1i * freq / 299792458) * ... 
%				(rodata.posITRF_fl * srcposITRF (upcals, :).'));
				(rodata.posITRF_fl * srcposhat (upcals, :).'));

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
	fit.peakfl = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * acc (:));
	fit.thsrc_wsf =  thsrc_wsf;
	fit.phisrc_wsf = phisrc_wsf;
	% fit.thsrc_wsf = racal; 
	% fit.phisrc_wsf = deccal;
