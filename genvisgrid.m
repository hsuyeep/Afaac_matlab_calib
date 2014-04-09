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
%		 .lim : Limit of convoutional kernel, in m.
%        .pa(1..4): Array of other parameters. For Gaussian GCF, 
%				  pa(1)/pa(2) = sigx/sigy of GCF gaussian.
%	deb   : Debug, turns on intermediate plotting etc.

function [gridvis] = genvisgrid (acc, u, v, parm, deb)
	
	% Generate UV grid
	gridrng = [-floor(parm.Nuv/2)*parm.duv:parm.duv:floor(parm.Nuv/2)*parm.duv];
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

	% plot statistics
	if (deb > 0)
		figure; 
		subplot (121);
		plot (uc, vc, gridvis, '.');
		subplot (122);
		plot (uc, vc, gridviscnt, '.');
	end;
