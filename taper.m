% Script to produce a visibility taper based on the parameters specified.
% The function incorporates the taper on the low spatial frequencies (short
% baseines) for spatial filtering, as well as on the regular uv plane.
% pep/20Aug13

% Arguments:
%   posITRF: Matrix of antenna positions in the array.
%   parm   : Structure containing parameters of the type of taper to apply.
%           type : ['Gaussian', 'Blackman']
%		  	pa(1/2): sigx/sigy of inner taper. 
%		  	pa(3/4): sigx/sigy of outer taper. 
%		  	minlambda: lower limit of visibility cut off(lambda).
%         	maxmeters: higher limit of visiblity cut off(m).
%	densitybins: turn on weighting by the density of visibilities.
%			-1 => don't turn on density weighting
%			any other value => number of bins in uvdist. space
%	 freq : The frequency of observation in Hz.
%   debug : Turn on debug messages.
%  Returns:
%	intap : The inner taper, normalized to 1.
%  outtap : The outer taper, normalized to 1.
% density : The spatial density of visibilities, for every uvdist.
%    mask : The final mask.
%  uvdist : The actual uvdistance of each visibility. 
% NOTE: If intap and outtap are being used, be sure to normalize the peak of the
% product to 1.

function [intap, outtap, density, mask, uvdist] = ...
								taper (posITRF, parm, densitybins, freq, debug)

	lambda = 299792458./freq; % In meters
	% NOTE: Assumes flagged antennas have been removed from posITRF.
	% ---- Generate uv coordinates (meters). ----
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	nants = length (posITRF); % Assumed symmetric

	u = meshgrid (posITRF(:, 1)) - meshgrid (posITRF(:, 1)).';
   	v = meshgrid (posITRF(:, 2)) - meshgrid (posITRF(:, 2)).';
    w = meshgrid (posITRF(:, 3)) - meshgrid (posITRF(:, 3)).';
    uvw = [u(:), v(:), w(:)];
	rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
    		   0.9928230000, -0.0954190000, 0.0720990000; ...
			   0.0000330000,  0.6030780000, 0.7976820000];
	uvw_rot = uvw * rotmat;   % All operations in the CS002 local plane
    % uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
    uvdist = sqrt(sum(uvw_rot.^2, 2) - (uvw_rot * normal).^2); % In meters

	% Generate density of baselines over uvdistance.
	density = ones (length (uvdist), 1);
	if (densitybins > 0)
		[n, bincen] = hist (uvdist, densitybins);
		for ind = 1:length (uvdist)
			[minval, index] = min(abs(bincen - uvdist(ind)));
			density(ind) = n(index);
		end;
		density = 1 - (density ./ max(density)); % Normalize
	end;

	%%%%%%%% test code for generating the inner ring taper on equispaced data.
	%%%%%%%% Change sigx/sigy for playing around with the taper.
	%	a = meshgrid ([-350:10:350]);
	%	b = a';
	%	ab = [a(:) b(:)];
	%	dist = sqrt (sum(ab.^2, 2));
	%	d = zeros (size (a));
	%	% Centered Gaussian
	%	d = exp (-(a.^2/(2*parm.pa(1)^2) + b.^2/(2*parm.pa(2)^2))); 
	%	figure;  subplot (121); imagesc (d); subplot (122); plot (dist, d(:));
	%
	%	% Sum of gaussians on the vertices of a square.
	%	cen = [-50 50];
	%	d = zeros (size (a));
	%	for ind = 1:2
	%		for jind = 1:2
	%			d = d + exp (-((a-cen(ind)).^2/(2*parm.pa(1)^2) +...
	%						   (b-cen(jind)).^2/(2*parm.pa(2)^2)));
	%		end;
	%	end;
	%	subplot (121); imagesc (d); subplot (122); plot (dist, d(:));
	%
	%	% Sum of gaussians on a circle of radius 50 m.
	%	d = zeros (size (a));
	%	circ = (sqrt (a.^2 + b.^2) == 50);
	%	circcoor = [a(circ) b(circ)];
	%	for ind = 1:length (circcoor)
	%		d = d + exp (-((a-circcoor(ind,1)).^2/(2*parm.pa(1)^2) + ...
	%					   (b-circcoor(ind,2)).^2/(2*parm.pa(2)^2)));
	%	end;
	%	subplot (121); imagesc (d); subplot (122); plot (dist, d(:));
	%
	%	% Surface of revolution of the function arctan(x), between 0 and 1.
	%	r = -pi:0.1:pi;
	%	z = atan (r);
	%	z = (z - min(z))./max (z-min(z)); % restrict to lie between 0 and 1.
	%	r = r + pi; % Move away from being symmetric around 0.
	%	theta = 0:pi/20:2*pi;
	%	xx = bsxfun(@times,r',cos(theta));
	%	yy = bsxfun(@times,r',sin(theta));
	%	zz = repmat(z',1,length(theta));
	%	surf(xx,yy,zz);
	%
	%	% Surface of revolution of a gaussian offset from the center.
	%	% x = -350:10:350;
	%	x = 150:10:350; % If the surface of revol. of half a gaussian is needed.
	%	z = exp (-(x-150).^2/(2*50^2));
	%	theta = 0:pi/20:2*pi;
	%	xx = bsxfun (@times, x', cos(theta));
	%	yy = bsxfun (@times, x', sin(theta));
	%	zz = repmat (z', 1, length(theta));
	%	surf (xx, yy, zz);
	
	    % Surface of revolution of a gaussian using another method for randomly 
		% sampled data.
	%	x = -350:10:350; 
	%	% y = -350:10:350; 
	%	y = linspace (-200, 200, length(x)); 
	%	sigx = 50; sigy = 50;
	%	[X, Y] = meshgrid (x, y);
	%	mask = zeros (size (X));
	%	xy = [X(:) Y(:)];
	%	dist = sqrt (sum(xy.^2, 2));
	%	distmat = reshape (dist, [length(x), length(x)]);
	%	xvec = X(:); yvec = Y(:);
	%	z = zeros (size (xvec));
	%	for ind = 1:length (xvec) % Do for all coord. pairs.
	%		% if ((dist(ind) > 50) && (dist(ind) < 200));
	%		if ((dist(ind) > 50))
	%			ang = atan2 (yvec(ind), xvec(ind)); % Angle of vector in radians
	%			an(ind) = ang;
	%			z(ind) = ...
	%
	%				exp(-( (xvec(ind) - 50*cos(ang)).^2/(2*50^2) +  ...
	%
	%					   (yvec(ind) - 50*sin(ang)).^2/(2*50^2)));
	%		end;
	%	end;
			
	%
	%%%%%%%% /test code for generating the inner ring taper on equispaced data.

	switch lower (parm.type)
		case {'gaussian'}
			fprintf (1, 'Taper mode: Gaussian, lambda: %f.\n', lambda);
			if (densitybins > 0)
				fprintf (1, 'Density weighting requested.\n');
			end;
			fprintf (1, ...
					'Inner sigx/sigy (lambda): %.2f/%.2f, outer:%.2f/%.2f.\n',...
					parm.pa(1), parm.pa(2), parm.pa(3), parm.pa(4));
			fprintf (1, 'At uv offset of %.2f/%.2f lambda/m.\n', ...
					 parm.minlambda, parm.minlambda*lambda);

			uvec = uvw_rot(:,1); vvec = uvw_rot(:,2);

			% Generate a 2D outer taper with specified params.
			% If variance < 0, implies no external taper is requested.
			if (parm.pa(3)  < 0 || parm.pa(4) < 0)
				outtap = ones (length (uvdist), 1);
			else
				outtap = exp (-(uvec.^2/(2*(parm.pa(3)*lambda)^2) + ...
								vvec.^2/(2*(parm.pa(4)*lambda)^2))); 
			end;

		 	% Generate an inner taper with specified parameters.
			% Gaussians on an idealized circle of radius parm.minlambda
			if ((parm.pa(1) > 0) && (parm.pa(2) > 0))
				intap = zeros (length (uvdist), 1);

				for ind = 1:length(uvec)
					if (uvdist(ind) > parm.minlambda*lambda)
						% Angle of vector in rads
						ang = atan2 (vvec(ind), uvec(ind)); 
						intap(ind) = ...
						exp(-( (uvec(ind)-lambda*parm.minlambda*cos(ang)).^2 /...
							   (2*parm.pa(1).^2) +  ...
							   (vvec(ind)-lambda*parm.minlambda*sin(ang)).^2 /...
							   (2*parm.pa(2).^2)));
					end;
				end;
				intap = 1-intap; 
				% set shortest baselines to 0;
				intap (uvdist < lambda*parm.minlambda) = 0; 
			else
				intap = ones (length (uvdist), 1);
			end;

		case {'arctan'}
			fprintf (2, 'Taper mode Arctan not implemented!\n');
		
		case {'triangle'}
			fprintf (2, 'Taper mode Triangle not implemented!\n');

		case {'hanning'}
			fprintf (2, 'Taper mode Hanning not implemented!\n');

		case {'hamming'}
			fprintf (2, 'Taper mode Hamming not implemented!\n');

		case {'blackman'}
			fprintf (2, 'Taper mode Blackman not implemented!\n');

		case {'notaper'}
			fprintf (2, 'No taper requested.\n');
			intap = ones (length (uvdist), 1);
			outtap = intap;
			density = intap;
		otherwise,
			fprintf (2, ...
			'Taper mode %s not known! Only applying short baseline filter.');
			tap = ones (nants);
	end;

	% -1 implies no inner taper
	%  0 implies hard inner taper at specified limits.
	% pep/21Apr14
	if ((parm.pa (1) == 0) || (parm.pa(2) == 0)) 
		fprintf ('taper: Applying hard spatial filter at min of %f lambda or % m.\n',...
				 parm.minlambda, parm.maxmeters);
		t = reshape (uvdist, [nants, nants]);
		% Hard cutoff mask.
 		intap = t(:) > min([lambda*parm.minlambda, parm.maxmeters]); 
	end;



	% Thin circle using sampled visibilities.
	% circ = (t > parm.minlambda & t < parm.minlambda+1); 
	% incirc = (t < parm.minlambda); % All visibilities within the circle.

	% Just use an inverted gaussian as an inner mask
	%	d = 5*exp (-((u.*incirc).^2/(2*parm.pa(1)^2) + ...
	%				 (v.*incirc).^2/(2*parm.pa(2)^2))); 
	%	d (d>1) = 0;
	%	d = 1-d;
	
	%	
	%	% Pick the baseline coordinates which fall on the circle.
	%	oncirc = find (circ(:) == 1);
	%	vecu = u(:); vecv = v(:);
	
	% For sum of 2D gaussians with cent. on a theoretical circle on the uv plane
	%	d = zeros (size (incirc));
	%	for ind = 1:length (th)
	%		fprintf (1, 'center: %.2f,%.2f.\n',parm.minlambda*cos(th(ind)), ...
	%					parm.minlambda*sin(th(ind)));
	%	e=exp(-((u.*incirc-parm.minlambda*cos(th(ind))).^2/(2*parm.pa(1)^2)...
	% 			  +(v.*incirc-parm.minlambda*sin(th(ind))).^2/(2*parm.pa(2)^2)));
	%		d = d + e;
	%		if (debug > 3)
	%			plot3 (u(incirc==1), v(incirc==1), e(incirc==1), '.');
	%			pause (0.1);
	%			hold on;
	%		end;
	%	end;

	% For sum of gaussians with centers on circ formed by measured uv coords
	%	for ind = 1:10:length(oncirc)
	%		fprintf (1, 'gauss with center: %.2f,%.2f.\n', ...
	%				 vecv(oncirc(ind)), vecv(oncirc(ind)));
	%		d = d + exp(-((u.*incirc-vecu(oncirc(ind))).^2/(2*parm.pa(1)^2)...
	% 			   +  (v.*incirc-vecv(oncirc(ind))).^2/(2*parm.pa(2)^2)));
	%		plot3 (u(incirc==1), v(incirc==1), d(incirc==1), '.');
	%	end;

	mask = intap .* outtap .* density;

	% Pull up mask values to lie between 0 and 1
	if (strcmp (lower(parm.type), 'notaper')== 0)
		mask = (mask - min(min(mask))) ./ (max(max(mask - min(min(mask)))));
	end;
	% d (incirc(:) == 0) = 1; % Set values outside the circle to 1.


	if (debug > 3)
		figure;
		subplot (131); plot (uvdist(:), intap(:), '.');
		title ('Inner taper'); xlabel ('uv dist. (m)'); ylabel ('taper amp');
		subplot (132); plot (uvdist(:), outtap(:), '.');
		title ('Outer taper'); xlabel ('uv dist. (m)'); ylabel ('taper amp');
		subplot (133); plot (uvdist(:), mask (:), '.');
		title ('Final mask'); xlabel ('uv dist. (m)'); ylabel ('taper amp');

		figure; 
		subplot (131); plot3 (uvw_rot(:,1), uvw_rot(:,2), intap(:), '.');
		title ('Inner taper'); xlabel ('u dist. (m)'); ylabel ('v dist. (m)');
		zlabel ('taper amp');
		subplot (132); plot3 (uvw_rot(:,1), uvw_rot(:,2), outtap(:), '.');
		title ('Outer taper'); xlabel ('u dist. (m)'); ylabel ('v dist. (m)');
		zlabel ('taper amp');
		subplot (133); plot3 (uvw_rot(:,1), uvw_rot(:,2), mask (:), '.');
		title ('Final mask'); xlabel ('u dist. (m)'); ylabel ('v dist. (m)');
		zlabel ('taper amp');
		% hold on;
		% plot3 (u(circ==1), v(circ==1), d(circ==1), '.r');
		% title ('2D taper as sum of gauss, on approx circle');
	end;

