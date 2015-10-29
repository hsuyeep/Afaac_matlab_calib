% Script to generate the density of UV plane sampling of a particular array configuration
% pep/01May15
% Arguments:
%  posITRF : ITRF locations of the antennas of this array
%  gridsize: The gridsize in meters
%  weight  : The weight associated with each visibility.

function [gridvis] = genvisden (posITRF, gridsize, weight)

	% Generate visibilities
	u = meshgrid (posITRF(:,1)) - meshgrid (posITRF(:,1)).';
	v = meshgrid (posITRF(:,2)) - meshgrid (posITRF(:,2)).';
	w = meshgrid (posITRF(:,3)) - meshgrid (posITRF(:,3)).';
	u = u(:); v = v(:); w = w(:);
    
    assert (length(u) == length(weight));
    uvw = [u(:), v(:), w(:)];
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
    		   0.9928230000, -0.0954190000, 0.0720990000; ...
			   0.0000330000,  0.6030780000, 0.7976820000];
	uvw_rot = uvw * rotmat;   % All operations in the CS002 local plane
    uvdist = sqrt(sum(uvw_rot.^2, 2) - (uvw_rot * normal).^2); % In meters

	% Generate a grid on the UV plane. Using 2-D plane, as W component will be projected onto this anyway.
	range = linspace (-max(uvdist), max(uvdist), 2*max(uvdist)/gridsize);
	maxu = max(u);
	maxv = max(v);
	gridvis = NaN (ceil(2*maxu/gridsize), ceil(2*maxv/gridsize));
    nrows = size (gridvis, 1);
    ncols = size (gridvis, 2);
	% For every visibility, find the cell in which it falls, increment the cell contents by the appropriate weight
	for ind = 1:length(u)
		uind = floor(u(ind)/gridsize)+floor(nrows/2)+1;
		vind = floor(v(ind)/gridsize)+floor(ncols/2)+1;

        if (isnan(gridvis(uind, vind)))
            gridvis (uind, vind) = 1;
        else
            gridvis (uind, vind) = gridvis (uind, vind) + weight(ind);
        end;
	end;
