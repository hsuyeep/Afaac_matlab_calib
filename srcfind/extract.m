% Script to extract the fluxes and positions of sources specified from the 
% 3CR catalog, given a calibrated AARTFAAC image.
% NOTE: If the Sun is not subtracted, a point src. model can be subtracted.
% pep/0Jan13

function [fitparams, goodfit] = extract (calmap,l,m, tobs, freq, abovejy, catalog, ...
								sig, det_thresh, debug)
	tobs_jd = tobs/86400. + 2400000.5; % Convert to JD.
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	det_limit = det_thresh*sig;
	goodfit = 0;
	persistent maphdl;
	persistent modelhdl;

	% radius of region to consider for source extraction, in pixel units.
	rad = 5;	% A square of side 20 pixels.
	nrow = length(l); ncol = length(m);
	dl = (max(l) - min(l))/nrow;
	dm = (max(m) - min(m))/ncol;


	% longitude and latitude of CS002 (AARTFAAC center).
	lon = 6.869837540;                         % longitude of CS002 in degrees
	lat = 52.915122495;                        % latitude of CS002 in degrees 

	% NOTE: There is no overlap between A-team sources and LOFAR calibrators.
	Ateamsrcind = [324, 283, 88, 179];
	lofarcalsrcind = [34, 89, 117, 192, 200, 256];
	% Indices in 3CR catalog: 3c[48=34,147=89,196=117,286=192,295=200,380=70]
	lofarcalsrcs = catalog(lofarcalsrcind);
	Ateamsrcs = catalog(Ateamsrcind); 

	% Treat the sources chosen on the basis of a flux cutoff.
	sel = ([catalog(:).flux] > abovejy);
	selsecind = find (sel > 0);
	finalselind = union ([Ateamsrcind, lofarcalsrcind], selsecind);
	srclist = catalog (finalselind);
	% srclist = catalog (selsecind);

	% Determine how many of our sources are above the horizon
	srcpos = radectoITRF([srclist.alpha], [srclist.delta], true, tobs_jd);
	up = srcpos * normal > 0;
	srclist = srclist(up);
	% selsecind = selsecind (up); % Also update the selection of sources.
	finalselind = finalselind (up); % Also update the selection of sources.

	% Check for residues of A-team sources, eliminate if too weak a residual
	% is found.
	% Indices of Ateam sources in srclist.
	%	Ateamind = find (ismember (srclist, Ateamsrcs)); 
		
	%	for ind = 1:length(Ateamind)
	%		% Find a region around the A-team source
	%		[srcl, srcm] = radectolm(srclist(Ateamind(ind)).alpha, ... 
	%			srclist(Ateamind(ind)).delta, tobs_jd, lon, lat);
	%
	%		% NOTE:The 'l' or 'East-West' axis is along Y axis, or over rows,
	%		% NOTE:The 'm' or 'North-South' axis is along Xaxis, or over cols.
	%		srcy = round (srcl/dl) + nrow/2 + 1; 
	%		srcx = round (srcm/dm) + ncol/2 + 1;
	%
	%		
	%		% Compare the maximum pixel value with the global image rms
	%
	%		% If not satisfied, delete the Ateam source from the srclist.
	%	end;



	% Locate the approximate positions of sources of interest, based on the
	% catalog positions.
	[srcl, srcm] = radectolm([srclist.alpha], [srclist.delta], tobs_jd, ... 
								lon, lat);
	% NOTE: The 'l' or 'East-West' axis is along the Y axis, or over rows,
	% NOTE: The 'm' or 'North-South' axis is along the X axis, or over cols.
	srcy = round (srcl/dl) + nrow/2 + 1; 
	srcx = round (srcm/dm) + ncol/2 + 1;

	% Create images only if debug level is high enough
	if debug > 0 
		if isempty(maphdl) % No figure created for debug!
			maphdl = figure;
			modelhdl = figure;
		end
		figure (maphdl); imagesc (calmap);

		% Overplot all selected srcs
		for ind = 1:length(srclist)
		    text ('Color' , [1 1 1], 'Position', [srcx(ind) srcy(ind)], ... 
				  'String', srclist(ind).name);
		end;
	end;

	% NOTE: 5 parameters are being estimated + srclist index + sigma.
	srcfitparams = struct ('xpos',0, 'ypos',0, 'amp',0, 'sigx',0, ... 
							'sigy',0, 'thresh',0, 'srcind',0, 'mjd',0);
	fitparams = zeros (length (srclist), 7);

	% Carry out source extraction of selected sources
	for ind = 1:length(srclist)
		% FIXME: For now, just ignore A-team sources
		if (ismember(srclist(ind).name, {Ateamsrcs.name}))
			continue;
		end;

		% Check if source is within our bounding box.
		% TODO: Check for other conditions on bounding box.
		if srcx(ind)-rad < 1 | srcx(ind)+rad > nrow
			disp (['NOTE: Source ' srclist(ind).name ' too close to border']);
			continue;
		end;

		if srcy(ind)-rad < 1 | srcy(ind)+rad > nrow
			disp (['NOTE: Source ' srclist(ind).name ' too close to border']);
			continue;
		end;
		
		% Extract out the subregion corresponding to this src.
		srcyreg = [srcy(ind)-1*rad:srcy(ind)+1*rad];
		srcxreg = [srcx(ind)-1*rad:srcx(ind)+1*rad]; 
		if debug > 1
			disp ( ...
		sprintf ('Source name: %s, subimage range:X (%3d, %3d), Y (%3d, %3d)', srclist(ind).name, min(srcxreg), max(srcxreg), min(srcyreg), ...
			max(srcyreg)));
		end
		% NOTE: To extract a subarray, we need array(row_range, col_range),
		% NOTE: Since the Image Y-axis goes over rows, and X-axis over cols,
		% NOTE: we need srcyreg in the rows position.
		data = calmap (srcyreg, srcxreg); 

		% NOTE: Origin of images in pixel units is always top left corner.
		% NOTE: ROWS of an image correspond to the vertical Y-axis,
		% NOTE: COLS of an image correspond to the horizontal X-axis.
		% Find the position of the source in the whole image, in pixel units.
		% Y = row vec. holding max. across all rows of everycol., 
		% I = index of row holding largest value among rows of every column.
		[Y, I] = max (data); 

		% Gives the peak flux across all cols and rows of subimage, 
		% row of peak pixel in Icol (scalar).
		[peak_flux, peak_col] = max(Y); 
		peak_row = I(peak_col);
		peak_x = srcxreg(peak_col); % Now in pixel units.
		peak_y = srcyreg(peak_row); % Now in pixel units.
		
		if debug > 1
			figure; imagesc (srcxreg, srcyreg, data);
	    	text ('Color' , [1 1 1], 'Position', [srcx(ind) srcy(ind)], ... 
				  'String', srclist(ind).name);
		end;

		disp (sprintf ('Source name: %s. SNR: %f', srclist(ind).name,...
			 peak_flux/sig));

		% Also ignore sources which do not have a good SNR.
		if peak_flux < det_limit
			disp (sprintf ('<---Ignoring src: %s,peak flux: %f,limit: %f', ...
					srclist(ind).name, peak_flux, det_limit));
			continue;
		end;

		% Fit gaussian to subimage. Parameters passed are
		sigx = 1; sigy = 1;
		init_par = [peak_x, peak_y, peak_flux, sigx, sigy];
        goodfit = goodfit + 1;
		fitparams(goodfit,1:5) = ... 
				fit2dgauss(data,srcxreg,srcyreg,init_par,debug);
		fitparams(goodfit,6) = sig;
		% fitparams(goodfit,7) = selsecind(ind); % Index into the passed catalog.
		fitparams(goodfit,7) = finalselind(ind); % Index into the passed catalog.
		fitparams(goodfit,8) = tobs;
	
		
		disp (sprintf ('Catalog flux: %8.4f Est. peak: %8.4f, (x,y): (%03d, %03d), sig(x,y): (%6.4f, %6.4f)', srclist(ind).flux, fitparams(goodfit, 3), fitparams(goodfit, 1), fitparams(goodfit, 2), fitparams(goodfit, 4), fitparams(goodfit, 5))); 
		% pause;

	end 

	disp (sprintf ('Extracted %d models successfully out of %d in list', ... 
		  goodfit, length(srclist)));

	% Generate a model sky based on estimated parameters, also the difference
	% between the original and the model sky.
	if debug > 0
		model_sky = zeros (size (calmap));
		[lmesh, mmesh] = meshgrid ([1:size(calmap,1)], [1:size(calmap,2)]); 
		for ind = 1:goodfit
		 	model_sky = model_sky + fitparams(ind,3)* ...
				exp(-(lmesh-fitparams(ind,1)).^2/(2*fitparams(ind,4)^2)).* ...
				exp(-(mmesh-fitparams(ind,2)).^2/(2*fitparams(ind,5)^2));
		end
		figure (modelhdl);
		imagesc (model_sky);
	
		% Overplot model source names
		for ind = 1:goodfit
            disp (sprintf ('ind: %d, fitparam(ind,7): %d.', ind, fitparams(ind,7)));
		   	text ('Color', [1 1 1], 'Position', [fitparams(ind,1), ...
		fitparams(ind,2)], 'String', catalog(fitparams(ind,7)).name);
		end;	
	end
