% Script to carry out convolutional gridding of recorded visibilities
% pep/07Apr14
% Arguments:
%	 acc  : observed visibilities
%    u,v  : UV coordinates of observed visibilities, in m.
%    parm : Parameters associated with the convolutional kernel.
%        .type: Type of GCF; Currently only Gaussian supported
%   	 .duv : Grid spacing of GCF, in WAVELENGTH units. 
%    	 .Nuv : Grid extent, in number of grid cells.
%		 .uvpad: Padding to be added to UV grid, for higher resolution.
%		 .lim : Limit of convoutional kernel, in WAVELENGTH units.
%        .pa(1..4): Array of other parameters. For Gaussian GCF, 
%				  pa(1)/pa(2) = sigx/sigy of GCF gaussian, in WAVELENGTH units.
%   freq  : Frequency of observation in Hz.
%	deb   : Debug, turns on intermediate plotting etc.

%  Returns:
%  gridvis: The padded gridded visibility values of dim parm.uvpad x parm.uvpad
%  gridviscnt: The count of visibilities contributing to a single grid point.
%  padgridrng: The coordinates of the padded gridded visibilities, assumed symmetrical.
%  kern   : The applied FFT kernel.

function [gridvis, gridviscnt, padgridrng, kern] = genvisgrid (acc, u, v, parm, freq, deb)
% function [gridvis] = genvisgrid (acc, u, v, parm, freq, deb)
	
	if (isempty (acc) || isempty (parm))
		error ('genvisgrid: Cannot grid with no acc OR no parameters specified!');
	end;

	% Convert u,v coordinates of visibilities to wavelength units.
	lambda = 299792458./freq;
	u = u(:) ./ lambda; v = v(:) ./ lambda;

	% Generate UV grid in wavelength units
	uvlim = floor (parm.Nuv/2)*parm.duv;
	gridrng = linspace (-uvlim, uvlim, parm.Nuv);
	paduvlim = floor (parm.uvpad/2)*parm.duv;
	padgridrng = linspace (-paduvlim, paduvlim, parm.uvpad);

	%gridrng=[-floor(parm.Nuv/2)*parm.duv:parm.duv:(floor(parm.Nuv/2)-1)*parm.duv];
	[uc, vc] = meshgrid (gridrng);
	uc = uc (:); vc = vc (:); % Vectorize the grid.
	% ngrid = length (gridrng);

	% Generate convolutional kernel with the same gridsize as underlying grid.
	kernrng = [-parm.lim:parm.duv:parm.lim];
	[uk, vk] = meshgrid (kernrng);	
	fprintf (1, 'NOTE: GCF kernel size is %dx%d cells.\n', ...
			 length (kernrng), length(kernrng));

	% Currently adding visibility gridding in an ad-hoc manner
	switch (lower (parm.type))
		case 'gaussian'

		fprintf (1, 'Gaussian GCF requested.\n');
		if (deb > 0)
			fprintf (1, 'parm.duv=%.2f, parm.Nuv=%d, parm.lim=%.2f, parm.pa(1)=%.2f\n',...
					 parm.duv, parm.Nuv, parm.lim, parm.pa(1));
		end;

	
	
		% Initialize matrix for holding gridded visibilities and the intermediate
		% weights.
		gridvis = zeros (1, length (uc));
		gridviscnt = gridvis; 
		weight = zeros (size (acc));

		% Generate gcf (for illustration purposes only)
		kern = 1./(2*pi*parm.pa(1)*parm.pa(2)) * exp (-((uk.^2 / (2*parm.pa(1)^2)) + (vk.^2 / (2*parm.pa(2)^2))));
		% kern_scale_fact = sum(sum(kern));
		% kern = kern./kern_scale_fact; % Make sure kernel integral is 1
		fprintf (1, '<-- Convolutional kernel integral: %f\n', sum(sum(kern)));
	
		% Convolve visibilities with GCF 
		% For each grid point:
		for grvis = 1:length (uc(:))
			% Generate a list of visibilities falling within 
			% the GCF range centered at this grid point.
			vislist = find (sqrt ((u-uc(grvis)).^2 + (v-vc(grvis)).^2) < 2*parm.lim);
	
			% Generate the weights associated with these visibilities for the given
			% GCF theoretical kernel.
			for sel = 1:length(vislist)
				weight(vislist(sel)) ...
				 = (1./(2*pi*parm.pa(1)*parm.pa(2)))*exp (-(((u(vislist(sel)) - uc (grvis)).^2 / (2*parm.pa(1)^2)) +...
					      (v(vislist(sel)) - vc (grvis)).^2 / (2*parm.pa(2)^2)));
			end;
	
			% Store the weighted sum of these visibilities at the chosen grid point
			gridvis (grvis) = sum (acc (vislist) .* weight(vislist));
			weight (vislist) = 0;
	
			% Generate statistics
			% How many vis. contributed to this grid point?
			gridviscnt (grvis) = length (vislist); 
			if (mod (grvis, 1000) == 999) fprintf (1, '.'); end;
		end;

		% Pad gridded visibilities for higher resolution
	    % zero padding to desired (u,v)-size
	    N = parm.Nuv;
	    N1 = floor(parm.uvpad/2) - floor (N/ 2);
	    % N2 = ceil((parm.uvpad + 1 - N) / 2) - 1; % Removed, works for odd Nuv. pep/23Apr14
	    % N2 = ceil((parm.uvpad + 1 - N) / 2);
		N2 = parm.uvpad - (N1 + parm.Nuv);
	    
	    % Surround gridded visibilities with 0-padding to create padded visibility 
	    % matrix.
	    gridvis = [zeros(N1, parm.uvpad); ...
	               zeros(N, N1), reshape(gridvis, [parm.Nuv parm.Nuv]) zeros(N, N2); ...
	               zeros(N2, parm.uvpad)];
	    gridvis(~isfinite(gridvis)) = 0;
		% Vectorize again
		gridvis = gridvis (:);

	case 'pillbox'
		fprintf (1, 'Pillbox GCF requested.\n');
		parm.pa(1) = 0; parm.pa(2) = 0; 
		if (deb > 0)
			fprintf (1, 'parm.duv=%.2f, parm.Nuv=%d, parm.uvpad=%d, parm.lim=%.2f, parm.pa(1)=%.2f\n',...
					 parm.duv, parm.Nuv, parm.uvpad, parm.lim, parm.pa(1));
		end;

		% Generate gcf (for illustration purposes only)
		kern = ones (size (uk));

	    % create object for interpolation
	    gridvis = zeros(parm.Nuv);
		gridviscnt = gridvis;
		missed_vis = 0;   % cumulative count of ignored visibilities due to 
						  % gridded value exeeding grid size.
	    %W = zeros(Nuv);

        % Vectorised version of the pillbox gridding process
	    % Get amp. and direction vector of visibility
        ampl = abs (acc);
        phasor = acc ./ ampl;

	   	% Determine the grid along U-axis in which observed visibility falls.
        uidx = u ./ parm.duv + parm.Nuv/2;
        uidxl = floor (uidx);
        uidxh = ceil (uidx);

        % Find absolute distance of measured visibility from grid points, in the
    	% U-direction only.
	    dul = abs(uidx - uidxl);		
	    duh = abs(uidx - uidxh);

        % Distribute the visiblity amplitude among the two grid points on the 
		% U-axis in proportion to their distance from the observed visiblity.
	    sul = duh .* ampl;
	    suh = dul .* ampl;
        
    	% Determine the grid along V-axis in which observed visibility falls.
	    vidx = v ./ parm.duv + parm.Nuv / 2;
	    vidxl = floor(vidx);		% Find the lower and higher gridded V-value.
	    vidxh = ceil(vidx);

        % Find absolute distance of measured visibility from grid points, in the
    	% V-direction only.
	    dvl = abs(vidx - vidxl);
	    dvh = abs(vidx - vidxh);

    	% Distribute the HIGHER u-grid point's share of the observed visibility
		% amp. between the higher and lower V-grid point.
	    sull = dvh .* sul;
	    suhl = dvh .* suh;
	    
	    % Distribute the LOWER u-grid point's share of the observed visibility 
		% amp. between the higher and lower V-grid point.
	    sulh = dvl .* sul;
	    suhh = dvl .* suh;

    	% Now that the observed visiblity amplitude is distributed among its 
    	% surrounding 4 grid points, fill the gridded visibility matrix with
    	% vectors with the same phase as the original observed visibility.
    	% NOTE: Adding the 4 vectors at the corners of the grid square will give
    	% back the original ungridded observed visibility.
    	% NOTE: We need to accumulate this to prevent overwriting the gridded 
		% values from a visibility from a neighbouring grid square.
        missed_vis =  sum (uidxl < 1)        + sum (uidxh < 1)       + ...
                      sum (uidxl > parm.Nuv) + sum (uidxh > parm.Nuv)+ ...
                      sum (uidxl < 1)        + sum (uidxh < 1)       + ...
                      sum (uidxl > parm.Nuv) + sum (uidxh > parm.Nuv);


		% Deal with the autocorrelations explicitly, else they are set to 0 and lost.
        ac_ind = (u == 0 && v == 0)
		gridvis(ac_ind) = gridvis(ac_ind) + ampl (ac_ind);

        lin_ind_ll = sub2ind (size (gridvis), uidxl, vidxl);
        lin_ind_lh = sub2ind (size (gridvis), uidxl, vidxh);
        lin_ind_hl = sub2ind (size (gridvis), uidxh, vidxl);
        lin_ind_hh = sub2ind (size (gridvis), uidxh, vidxh);

        gridvis (lin_ind_ll) = gridvis (lin_ind_ll) + sull .* phasor;
        gridvis (lin_ind_lh) = gridvis (lin_ind_lh) + sulh .* phasor;
        gridvis (lin_ind_hl) = gridvis (lin_ind_hl) + suhl .* phasor;
        gridvis (lin_ind_hh) = gridvis (lin_ind_hh) + suhh .* phasor;

        gridviscnt (lin_ind_ll) = gridviscnt (lin_ind_ll) + 1;
        gridviscnt (lin_ind_lh) = gridviscnt (lin_ind_lh) + 1;
        gridviscnt (lin_ind_hl) = gridviscnt (lin_ind_hl) + 1;
        gridviscnt (lin_ind_hh) = gridviscnt (lin_ind_hh) + 1;

		if (missed_vis > 0)
		 	fprintf (2, 'Missed vis: %d\n', missed_vis); 
		end;
	
	    % zero padding to desired (u,v)-size
	    N = parm.Nuv; % size (gridvis, 1);
	    N1 = floor(parm.uvpad/2) - floor (N/ 2);
		N2 = parm.uvpad - (N1 + parm.Nuv);
	    % N1 = floor((parm.uvpad - N) / 2);
	    % N2 = ceil((parm.uvpad + 1 - N) / 2) - 1;
	    
	    % Surround gridded visibilities with 0-padding to create padded visibility 
	    % matrix.
	    gridvis = [zeros(N1, parm.uvpad); ...
	               zeros(N, N1), gridvis, zeros(N, N2); ...
	               zeros(N2, parm.uvpad)];
	    gridvis(~isfinite(gridvis)) = 0;

		% Generate a gridded uv coordinate for display purposes
		% gridrng = [-floor(parm.uvpad/2)*parm.duv:parm.duv:(floor(parm.uvpad/2)-1)*parm.duv];
		% [uc, vc] = meshgrid (gridrng);
		% uc = uc (:); vc = vc (:); % Vectorize the grid.
		gridvis = gridvis (:);  % Vectorize padded and gridded visibilities

	case 'exponential'
		error ('genvisgrid: GCF Exponential not yet implemented!');

	case 'sinc'
		error ('genvisgrid: GCF Sinc not yet implemented!');

	case 'expsinc'
		error ('genvisgrid: GCF Exponential*Sinc not yet implemented!');

	case 'spheroid'
		error ('genvisgrid: GCF Spheriod not yet implemented!');

	otherwise
		error (2, 'genvisgrid: Unknown gridding type! Quitting!\n');
		
	end;

	% plot statistics
	if (deb > 0)
		figure; 
		subplot (221);
		plot3 (u, v, abs(acc), '.');
		colorbar;
		view (0, 90);
		xlim ([gridrng(1) gridrng(end)]); % Assumed symmetrical grid.
		ylim ([gridrng(1) gridrng(end)]);
		xlabel ('u(\lambda)'); ylabel ('v(\lambda)'); title ('Ungridded visibilities');

		subplot (222);
		imagesc (padgridrng, padgridrng, reshape (abs(gridvis), [parm.uvpad parm.uvpad])); % Assumed symmetrical grid.
		colorbar;
		xlabel ('u(\lambda)'); ylabel ('v(\lambda)'); title ('Gridded visibilities');

		
%		plot3 (uc, vc, abs(gridvis(:)), '.');
%		xlabel ('u (m)'); ylabel ('v(m)'); zlabel ('Gridded vis. value');
%		subplot (122);
%		plot3 (uc, vc, gridviscnt(:), '.');
%		xlabel ('u (m)'); ylabel ('v(m)'); zlabel ('No. of contributing vis. per grid point');
	
		% Plot the GCF and its FFT
		% Generate the GCF
		subplot (223);
		gcf = zeros (length (padgridrng));
		ind = int32 (parm.uvpad/2 - length (kernrng)/2);
		gcf (ind:ind+length(kernrng)-1, ind:ind+length(kernrng)-1) = kern;
		imagesc (gridrng, gridrng, gcf); colorbar;
		xlabel ('u(\lambda)'); ylabel ('v(\lambda)'); 
		title ('GCF, visibility domain');
    	gcf = conj(flipud(fliplr(fftshift(gcf))));
		
		subplot (224);
		dl = 1/(parm.uvpad*parm.duv);
		lmax = dl * parm.uvpad/2;
		lrng = linspace (-lmax, lmax, parm.uvpad);
    	GCF = fftshift(fft2(gcf));
		GCF = GCF ./ (max(max(abs(GCF)))); % Is this normalization valid??

		% Restrict image to lie between -1 and 1 l,m coordinates.
		d = intersect (find (lrng > -0.2), find (lrng < 0.2)); 
		% imagesc (lrng(d), lrng(d), abs(GCF(d, d))); colorbar;
		plot (lrng (d), abs(GCF (int32(length (GCF)/2), d)));
		
		xlabel ('l'); ylabel ('m'); title ('GCF, Image domain');
		
		figure;
    	gridvis = conj(flipud(fliplr(fftshift(reshape (gridvis, [parm.uvpad parm.uvpad])))));
    	psf = fftshift(fft2(gridvis));
		% psf = fftshift(fft2(reshape(gridvis, [parm.uvpad parm.uvpad])));
		psf = psf ./ max(max(abs(psf)));

		% Carry out grid correction on the generated PSF.
		psf = psf./GCF;

		subplot (221);
		plot (lrng(d), 20*log10(abs(psf (length(psf)/2, d))));
		xlabel ('l'); ylabel ('Power (dB)');
		ylim ([-75 0])
		% xlim ([lrng(1) lrng(end)]);
		grid on;
		subplot (223);
		plot (lrng(d), 20*log10(abs(psf (d,length(psf)/2))));
		xlabel ('m'); ylabel ('Power (dB)');
		ylim ([-75 0])

		% xlim ([lrng(1) lrng(end)]);
		grid on;
		subplot (122);
		imagesc (lrng(d), lrng(d), 20*log10(abs(psf(d,d))));
		xlabel ('l'); ylabel ('m'); title ('FFT image');
		colorbar;
		caxis ([-75 0]);
	end;

	% Reshape gridded visibility to maintain compatibility with fft_imager*
	gridvis = reshape (gridvis, [parm.uvpad parm.uvpad]);
