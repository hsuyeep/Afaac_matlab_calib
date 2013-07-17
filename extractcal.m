% Script to extract out the peak and integrated flux of calibrator sources seen 
% in a given snapshot image.
% pep/12Jul13
%  Arguments: 
%	 img  : Image on which to operate
%     l,m : l,m coordinates to go with image.
%	 tobs : Timestamp of snapshot in MJD seconds
%    freq : Frequency of observation in Hz.
% catalog : Specified catalog containing positional information of cal sources.
% callist : Calibrator list as indices within catalog.
%   epoch : True if positions are B1950, false if they are J2000.
%     psf : Matrix containing the estimated PSF for the current image.
%			If empty, a 2D gaussian is fit to sources.
%   debug : If 1, display island maps, fit residuals/curves, more msgs.

% Returns:
%  upcals : Bool vector containing 1s for the calibrators actually used fron
%			specified callist. 
%    fit  : An array of structures of size sum (upcals), containing the fit 
%			parameters as returned by the underlying fitting routine.
% Members of fit:
% fitparam: Vector containing optimization results of fit parameters. Check
%			underlying fitting routine for definition of members.	  
%   resid : The norm of the residual between model and data. Use this to 
%			determine fit quality.
%  exitfl : The exit status of the func. minimization. 0 implies things are good

function [upcals, fit] = ...
		extractcal (img, l, m, tobs, freq, rodata, catalog, callist, epoch,...
					 psf, debug)


	% addpath ../
	% addpath srcfind/;

    tobs_jd = tobs/86400 + 2400000.5; % Convert from MJD secs. to JD day units
	errrad = 5;  % In pixel units.
	fitrad = 5; % In pixel units

	% NOTE: This assumes image is square.
	N = size (img, 1); % Number of pixels in image.
	lm2pix  = (max(l) - min(l))/N;


	% Extract out calibrators within FoV and above horizon
	racal  = [catalog(callist).alpha];
	deccal = [catalog(callist).delta];
	fluxcal= [catalog(callist).flux];
	
	srcposITRF = radectoITRF (racal, deccal, epoch, tobs_jd);
	upcals = srcposITRF * rodata.normal > 0;
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

	% Handle debug
	if debug > 2
		islemap = figure; 
		fit2dmap = figure;
	else
		islemap = -1;
		fit2dmap = -1;
	end;
	
	% For each calibrator
	for ind = 1:sum(upcals)	

		% Locate calibrator pixel island within error circle
		% Convert l,m to pixel coordinates
		% NOTE: l is East-west, but it goes along the rows, or the Y-axis of 
		% an imagesc image.
		rowpix = int32 (lup (ind)/lm2pix + N/2); 
		colpix = int32 (mup (ind)/lm2pix + N/2);
		
		% Search for peak in l,m within error square
		errsq = img (rowpix-errrad:rowpix+errrad, colpix-errrad:colpix+errrad);
		[colmax, colmaxind] = max (errsq);
		[rowmax, rowmaxind] = max (colmax);

		% The local location of the peak and its value.
		peakm = colmaxind (rowmaxind); peakl = rowmaxind;  % l = rows;
		peakpix = errsq(peakl, peakm);
		

		% Create pixel island of appropriate size 
		% From the position of observed max, choose a square of side/2=fitsquare
		% NOTE: Need to determine till where symmetry is maintained. The effect
		% of sidelobes can make a gaussian profile non-symmetric, causing 
		% problems in gaussian fitting.

		% convert max position to image coords.
		limg = rowpix - errrad + peakl - 1; % l is rows!
		mimg = colpix - errrad + peakm - 1; % m is cols!

		% Go along l-axis (so rows), checking for symmetry
		residth = 0.25; % In percent wrt. deviation by peak.
		resid = 0; lsig = 0;
		while resid < residth
			resid = abs (img(limg+lsig, mimg) - img(limg-lsig, mimg))/peakpix;
			lsig = lsig + 1;
		end;

		% Go along m-axis (so cols), checking for symmetry
		resid = 0; msig = 0;
		while resid < residth
			resid = abs (img(limg, mimg+msig) - img(limg, mimg-msig))/peakpix;
			msig = msig + 1;
		end;
			
		fitl = [limg-fitrad:limg+fitrad]; 
		fitm = [mimg-fitrad:mimg+fitrad];

		% NOTE: Currently not using this!
		% fitl = [limg-lsig:limg+lsig]; 
		% fitm = [mimg-msig:mimg+msig];
		lsig = 3; msig = 3;
		fitmat = img (fitl, fitm);
		fitl = double (fitl); fitm = double (fitm);

		% display offset from catalog l,m
		if debug > 1
			fprintf ('Catalog peak: (%d, %d), datapeak: (%d, %d), val: %f\n',...
				  	 rowpix, colpix, limg, mimg, errsq(peakm, peakl));
		end;

		% Fit model to calibrator island
		% init_par = [peak_x, peak_y, peak_flux, sigx, sigy];
		% NOTE: Conversion to double is required by the exp and fminsearch 
		% functions.
		init_par = double ([limg, mimg, errsq(peakl, peakm), lsig, msig]);

		if (isempty (psf))
			% NOTE: l is the x-axis, ie, columns!
			fit(ind) = fit2dgauss (fitmat, fitl, fitm, ...
									double(init_par), debug, fit2dmap);
		else
			fit(ind) = fit2dpsf (fitmat, fitl, fitm, ...
									double(init_par), debug, fit2dmap);
		end; 
	
		% if (debug > 0)
		%	fprintf (1, 'Src: %s, pk: %.2f, int: %.2f, Resid: %f.\n', ...
  		%			 cat(ind).name, peak (ind), int (ind), fitres(ind));
		% end; 

		if (debug > 2)
			figure (islemap);
			subplot (122);
			imagesc (limg, mimg, fitmat);
			title (sprintf ('Cal %d: %s, %d Jy, %.1f counts\n', ...
				ind, cat(ind).name, cat(ind).flux, errsq(peakl, peakm)));
			xlabel ('l'); ylabel ('m');

			subplot (221); 
			hold off;
			plot (fitl, max(fitmat ), 'b.-'); 
			hold on;
			plot (fitl, ...
			fit(ind).fitparams(3)*...
	exp(-(fitl-fit(ind).fitparams(1)).^2/(2*fit(ind).fitparams(4)^2)),'r.-');
			title ('max over cols'); xlabel ('l');

			subplot (223); 
			hold off;
			plot (fitm, max(fitmat'), 'b.-'); 
			hold on;
			plot (fitm, ...
			fit(ind).fitparams(3)*...
	exp(-(fitm-fit(ind).fitparams(2)).^2/(2*fit(ind).fitparams(4)^2)),'r.-');
			title ('max over rows'); xlabel ('m');
		end;
	end;
