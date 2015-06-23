% Script to generate a single, all-sky map from calibrated visibilities 
% per timeslice, available from a calsol file. Basically works by imaging
% the visibilities using a Direct Fourier Transform in the direction of the
% defined l,m grid.
% pep/17Mar15
% A modification of Stefan Wijnholds' cs10survey.m
% Arguments:
%    fname : File name of visibility file. Can be calibrated/uncalibrated.
%	 skip  : Number of records to skip.
%   offset : Initial offset in image file to start of processing.
%     nrec : Number of images to accumulate into all-sky map.
%			 -1 => operate on all images.
%      deb : Turn on debugging information.
%
%  Returns :
%		im0: dirty image;
%    skymap: Combined map;
%     alpha: RA coordinates of combined map.
%    delta : DEC coordinates of combined map.
%   the_res: The chosen resolution in elevation, for gridding the RADEC space.
%   phi_res: The chosen resolution in azimuth, for gridding the RADEC space.

function [im0,imcl, M, alpha, delta, the_res, phi_res] = allskymap(fname, skip, offset, nrec, flagant, deb)

	% Open visibility binary file.
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('allskymap: fid < 0! Quitting.');
		return;
	end;

	% Move by required offset
	for ind = 1:offset
    	try
			[acc, tobs, fobs] = readms2float(fid, -1, -1, 288);

		catch err
			fprintf (2, 'allskymap: EoF for visibility file reached!\n');
			return;
		end;
	end;

    try
		[acc, tobs, fobs] = readms2float(fid, -1, -1, 288);
	catch err
		fprintf (2, 'allskymap: EoF for visibility file reached!\n');
		return;
	end;
    
	if (skip == 0 || isempty (skip)) skip = 1; end;

	% Read in a record, moveby desired offset    
%	if (offset > 1)
%		fprintf (2, 'Moving to calib. solution offset %d\n', offset);		
%		for ind = 1:offset
%			rec = read(fid);
%		end;
%	end;

	if nrec < 0 
		[nrec, tmin, tmax, dt] = getnrecs (fname);
		fprintf (1, 'Found %d records at %f sec. resolution.', ... 
				nrec, dt);
	end;

	% Define parameters for gridding in RA/DEC space
	load ('poslocal_outer.mat', 'posITRF', 'poslocal');
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	uvdist = sqrt(uloc.^2 + vloc.^2);
	D = max(uvdist(:));
	c = 2.99792458e8;
	lambdamin = c/fobs;
	res = 0.75 * lambdamin / D;    % RA/DEC grid resolution, in rad.
    % res = 1.5 * lambdamin / D;    % RA/DEC grid resolution, in rad.	
	nlim = 0.025;
    cs002_lon = 6.869837540;  % longitude of CS002 in deg
    cs002_lat = 52.915122495; % latitude of CS002 in deg
		

    Rpos = [poslocal(:,1), poslocal(:,2), zeros(size(poslocal(:,1)))];

	% define a grid with resolution proportional with minimum projection
	limp = [];
	mimp = [];
	nimp = [];
	theta = cs002_lat*pi/180;
	phi_res = [];
	the_res = [];
	count = 1;
	while theta < pi/2
	    Nphi = ceil(2 * pi * cos(theta) / res);       
	    dphi = 2 * pi / Nphi;
		phi_res(count) = Nphi;
	    for phi_idx = 0:Nphi - 1
	        limp = [limp; cos(theta) * cos(phi_idx * dphi)];
	        mimp = [mimp; cos(theta) * sin(phi_idx * dphi)];
	        nimp = [nimp; sin(theta)];
	    end
	    theta = theta + res / max([abs(cos(theta - cs002_lat)), 0.2]);
		the_res(count) = res / max([abs(cos(theta - cs002_lat)), 0.2]);
		count = count + 1;
	end
	theta = cs002_lat*pi/180 - res;
	while theta > -pi/2
	    Nphi = ceil(2 * pi * cos(theta) / res);     
	    dphi = 2 * pi / Nphi;
		phi_res(count) = Nphi;
	    for phi_idx = 0:Nphi - 1
	        limp = [limp; cos(theta) * cos(phi_idx * dphi)];
	        mimp = [mimp; cos(theta) * sin(phi_idx * dphi)];
	        nimp = [nimp; sin(theta)];
	    end
	    theta = theta - res / max([abs(cos(theta - cs002_lat)), 0.2]);
		the_res(count) = res / max([abs(cos(theta - cs002_lat)), 0.2]);
		count = count + 1;
    end
	fprintf (1,'<-- Dec pixels: %d.\n', count);
    fprintf (1,'<-- Pixel dimensions of full map (l/m/n): %d %d %d\n', length(limp), length(mimp), length(nimp));
	Lgrid = [limp, mimp, nimp];
	alpha = atan2(limp, mimp) + pi;
	delta = asin(nimp);

	if (deb > 0)
		% [al_lin, de_lin] = meshgrid (linspace (-pi, pi, 1024), linspace (-pi/2, pi/2, 512));
   		[al_lin, de_lin] = meshgrid (linspace (0, 2*pi, 1024), linspace (-pi/2, pi/2, 512));        
	end;
	
	% The dirty image
	im0 = zeros (size (alpha));
	imcl = im0;

	% The count of the number of maps being accumulated in each pixel.
	nmap2pix = zeros (size (alpha));

% 	l = zeros(nrecs, size(Lgrid, 1));
% 	m = zeros(nrecs, size(Lgrid, 1));
 	M = 0;
	
	try
		for idx = offset:skip:offset+nrec
	        try
				[acc, tobs, fobs] = readms2float(fid, -1, -1, 288);
	        catch err
	             fprintf (2, 'EoF reached!');
	             rethrow(err);
	                    % break;
	        end;
	
			% Find flagged antennas
			badant = find (isnan (diag(acc)) | diag(acc) == 0);
			currflagant = union (badant, flagant);
			[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, currflagant); 
	    	antmask = zeros (size (acc));
	    	posmask = zeros (size (Rpos,1), 1);
	    	rem_ants = length(acc) - length(currflagant);
	    	for ind = 1:length(currflagant)
	    		antmask (currflagant(ind), :) = 1; antmask (:,currflagant(ind)) = 1;
	    	 	posmask (currflagant(ind), :) = 1;
	    	end
	    	acc = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
			
	        fprintf (1, '<-- rec %04d, %s, flagant: %s ', idx, datestr(mjdsec2datenum(tobs)), num2str(currflagant'));
	        tobs_jd = tobs/86400 + 2400000.5; % Convert MJDsec to JD.
			% determine the RA/dec grid points which fall within this obsevation.
			[l, m] = radectolm(alpha, delta, tobs_jd, cs002_lon, cs002_lat, 0);
	        up = (~isnan(l) & (sqrt(1 - l.^2 - m.^2) > nlim));
			inim = sum(up, 1) > 0;
	        l = l(up);
	        m = m(up);
	        n = sqrt(1 - l.^2 - m.^2);
	        Lgridsnapshot = [l, m, n];
	        insnapshot = (~isnan(l) & n > nlim);
			
	%         up = (~isnan(l) & (sqrt(1 - l.^2 - m.^2) > nlim));
	% 		inim = sum(up, 1) > 0;
	        
	        invR = inv(acc);
	        A = exp(-(2 * pi * i * fobs / c) * Rpos(posmask == 0,:) * Lgridsnapshot.');
	        A(:, ~insnapshot) = 0;
	        condR = cond(abs(invR).^2);
	        if condR > 1e10 || isnan(condR)
				fprintf (2, '<-- Condition number > 1e10!\n');
	            continue;
	        end
	        imdirty = (conj(A)' * conj(invR) .* A') * ones(sum(rem_ants), 1) - abs(A' * invR).^2 * inv(abs(invR).^2) * diag(invR);
			fprintf (1, 'Dirty pixels: %d.\n', length(imdirty));
	        Msnapshot = abs(A' * invR * A).^2 - abs(A' * invR).^2 * abs(invR).^2 * (invR * A).^2;
			% imclean = imdirty \ Msnapshot;
			imclean = Msnapshot \ imdirty;
	        im0(up == 1) = im0 (up == 1) + imdirty;
			imcl(up==1) = imcl(up == 1) + imclean';
	        % im0(up == 1) = im0 (up == 1) + imclean;
			nmap2pix (up == 1) = nmap2pix (up == 1) + 1;
	        % M(up == 1)  = M (up == 1) + Msnapshot;
	
			% Display the dirty map, if debugging is required.
			if (deb > 0)
				% radecmap = griddata (alpha, delta, im0, al_lin, de_lin, 'cubic');
                radecobj= TriScatteredInterp (alpha, delta, im0);
                radecmap = radecobj(al_lin, de_lin);
				imagesc (al_lin(1,:)*12/pi, de_lin(:,1)*180/pi, abs(radecmap)); colorbar;
                xlabel ('RA(Hr)'); ylabel ('Dec(deg)');
                title (sprintf ('Time: %s, Freq: %.2f\n', datestr(mjdsec2datenum(tobs)), fobs));
                drawnow();
			end;
	
	        % Carry out DFT imaging for each of the selected l,m. Accumulate
	        % the generated pixel value onto the RA/DEC grid.
	        % map = acm2skyimage (rec.calvis, posITRF(1,:), posITRF(2,:), rec.fobs, l(up == 1), m(up == 1));
	
			for ind = 1:skip
				[acc, img.tobs, img.freq] = readms2float (fid, -1, -1, 288);
			end;
	        
		end
	catch
		fprintf (1, '<-- End of input file reached.\n');
	end;

	im0 = im0 ./ nmap2pix;
    fclose(fid);
