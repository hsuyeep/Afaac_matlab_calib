% Script to carry out convolutional gridding of recorded visibilities
% pep/07Apr14
% Arguments:
%	 acc  : observed visibilities
%    u,v  : UV coordinates of observed visibilities, in m.
%    type : Type of convolutional kernel. 
%    parm : Parameters associated with the convolutional kernel.
%        .type: Type of GCF; Currently only Gaussian supported
%   	 .duv : Grid spacing of GCF, in m
%    	 .Nuv : Grid extent, in number of grid cells.
%		 .uvpad: Padding to be added to UV grid, for higher resolution.
%		 .lim : Limit of convoutional kernel, in m.
%        .pa(1..4): Array of other parameters. For Gaussian GCF, 
%				  pa(1)/pa(2) = sigx/sigy of GCF gaussian.
%	deb   : Debug, turns on intermediate plotting etc.

function [gridvis] = genvisgrid (acc, u, v, parm, deb)
	
	if (isempty (acc) || isempty (parm))
		error ('genvisgrid: Cannot grid with no acc OR no parameters specified!');
	end;

	% Currently adding visibility gridding in an ad-hoc manner
	switch (lower (parm.type))
		case 'gaussian'
		% Generate UV grid
		gridrng = [-floor(parm.Nuv/2)*parm.duv:parm.duv:(floor(parm.Nuv/2)-1)*parm.duv];
		[uc, vc] = meshgrid (gridrng);
		uc = uc (:); vc = vc (:); % Vectorize the grid.
		ngrid = length (gridrng);
	
		% Generate convolutional kernel
		kernrng = [-parm.lim:parm.duv:parm.lim];
	
		% Initialize matrix for holding gridded visibilities and the intermediate
		% weights.
		gridvis = zeros (1, length (uc));
		gridviscnt = gridvis; 
		weight = zeros (size (acc));
	
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
				 = exp (-((u(vislist(sel)) - uc (grvis)).^2 / (2*parm.pa(1)^2)) +...
					     ((v(vislist(sel)) - vc (grvis)).^2 / (2*parm.pa(2)^2)));
			end;
	
			% Store the weighted sum of these visibilities at the chosen grid point
			gridvis (grvis) = sum (acc (vislist) .* weight(vislist));
			weight (vislist) = 0;
	
			% Generate statistics
			% How many vis. contributed to this grid point?
			gridviscnt (grvis) = length (vislist); 
		end;

	case 'wijnholds'
	    % create object for interpolation
	    gridvis = zeros(parm.uvpad);
		gridviscnt = gridvis;
		missed_vis = 0;   % cumulative count of ignored visibilities due to 
						  % gridded value exeeding grid size.
	    %W = zeros(Nuv);
	    for idx = 1:length(u(:))            % For every recorded visibility
	    
	    	% Get amp. and direction vector of visibility
	        ampl = abs(acc(idx));			
	        phasor = acc(idx) / ampl;
	    
	    	% Determine the grid along U-axis in which observed visibility falls.
	        uidx = u(idx) / parm.duv + parm.Nuv / 2;  
	        uidxl = floor(uidx);		% Find the lower and higher gridded U-value.
	        uidxh = ceil(uidx);
	    
	        % Find absolute distance of measured visibility from grid points, in the
	    	% U-direction only.
	        dul = abs(uidx - uidxl);		
	        duh = abs(uidx - uidxh);
	    
	        % Distribute the visiblity amplitude among the two grid points on the 
			% U-axis in proportion to their distance from the observed visiblity.
	        sul = duh * ampl;
	        suh = dul * ampl;
	        
	    	% Determine the grid along V-axis in which observed visibility falls.
	        vidx = v(idx) / parm.duv + parm.Nuv / 2;
	        vidxl = floor(vidx);		% Find the lower and higher gridded V-value.
	        vidxh = ceil(vidx);
	    
	        % Find absolute distance of measured visibility from grid points, in the
	    	% V-direction only.
	        dvl = abs(vidx - vidxl);
	        dvh = abs(vidx - vidxh);
	    
	    	% Distribute the HIGHER u-grid point's share of the observed visibility
			% amp. between the higher and lower V-grid point.
	        sull = dvh * sul;
	        suhl = dvh * suh;
	    
	    	% Distribute the LOWER u-grid point's share of the observed visibility 
			% amp. between the higher and lower V-grid point.
	        sulh = dvl * sul;
	        suhh = dvl * suh;
	        
	    	% Now that the observed visiblity amplitude is distributed among its 
	    	% surrounding 4 grid points, fill the gridded visibility matrix with
	    	% vectors with the same phase as the original observed visibility.
	    	% NOTE: Adding the 4 vectors at the corners of the grid square will give
	    	% back the original ungridded observed visibility.
	    	% NOTE: We need to accumulate this to prevent overwriting the gridded 
			% values from a visibility from a neighbouring grid square.
			% fprintf (1, 'bline: %06d (%03.2f,%6.2f), uidx: %03d %03d, vidx: %03d, %03d\n', idx, u(idx), v(idx), uidxl, uidxh, vidxl, vidxh);
			if ((uidxl < 1) || (uidxh < 1) || (uidxl > parm.Nuv) || (uidxh > parm.Nuv))
				missed_vis = missed_vis + 1;
				continue;
			end;
	
			if ((vidxl < 1) || (vidxh < 1) || (vidxl > parm.Nuv) || (vidxh > parm.Nuv))
				missed_vis = missed_vis + 1;
				continue;
			end;
	
	        gridvis(uidxl, vidxl) = gridvis(uidxl, vidxl) + sull * phasor;
	        gridvis(uidxl, vidxh) = gridvis(uidxl, vidxh) + sulh * phasor;
	        gridvis(uidxh, vidxl) = gridvis(uidxh, vidxl) + suhl * phasor;
	        gridvis(uidxh, vidxh) = gridvis(uidxh, vidxh) + suhh * phasor;

			gridviscnt(uidxl, vidxl) = gridviscnt(uidxl, vidxl) + 1;
			gridviscnt(uidxl, vidxh) = gridviscnt(uidxl, vidxh) + 1;
			gridviscnt(uidxh, vidxl) = gridviscnt(uidxh, vidxl) + 1;
			gridviscnt(uidxh, vidxh) = gridviscnt(uidxh, vidxh) + 1;        
	        %W(uidx, vidx) = W(uidx, vidx)
   		 end
		if (missed_vis > 0)
		 	fprintf (2, 'Missed vis: %d\n', missed_vis); 
		end;
	
	    % zero padding to desired (u,v)-size
	    N = size(gridvis, 1);
	    N1 = floor((parm.uvpad - N) / 2);
	    N2 = ceil((parm.uvpad + 1 - N) / 2) - 1;
	    
	    % Surround gridded visibilities with 0-padding to create padded visibility 
	    % matrix.
	    gridvis = [zeros(N1, parm.uvpad); ...
	               zeros(N, N1), gridvis, zeros(N, N2); ...
	               zeros(N2, parm.uvpad)];
	    gridvis(~isfinite(gridvis)) = 0;

		% Generate a gridded uv coordinate for display purposes
		gridrng = [-floor(parm.uvpad/2)*parm.duv:parm.duv:(floor(parm.uvpad/2)-1)*parm.duv];
		[uc, vc] = meshgrid (gridrng);
		uc = uc (:); vc = vc (:); % Vectorize the grid.

	otherwise
		error (2, 'genvisgrid: Unknown gridding type! Quitting!\n');
		
	end;

	% plot statistics
	if (deb > 0)
		figure; 
		subplot (121);
		plot3 (uc, vc, abs(gridvis(:)), '.');
		subplot (122);
		plot3 (uc, vc, gridviscnt(:), '.');
	end;
