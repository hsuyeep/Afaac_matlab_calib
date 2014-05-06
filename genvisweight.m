% Script to generate the sweight for every visibility, based on the specified
% weighting scheme.
% pep/02Apr14
%  Arguments:
%   posITRF: The positions of the antennas making up the array
%   wparm  : The parameters to be passed onto the weighting function
%           wparm.type : 'Uniform', 'Superuniform (TODO)'
%			wparm.cellrad: The radius in meters in the uv space over which to 
%						   count for generating visibility weights.
%  plt : Bool to control generation of intermediate plots.
function [weight] = genvisweight (posITRF, wparm, plt)
	
	% Convert antenna positions into visibility values
	u = meshgrid (posITRF(:, 1)) - meshgrid (posITRF(:, 1)).';
   	v = meshgrid (posITRF(:, 2)) - meshgrid (posITRF(:, 2)).';
    w = meshgrid (posITRF(:, 3)) - meshgrid (posITRF(:, 3)).';
    uvw = [u(:), v(:), w(:)];
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
    		   0.9928230000, -0.0954190000, 0.0720990000; ...
			   0.0000330000,  0.6030780000, 0.7976820000];
	uvw_rot = uvw * rotmat;   % All operations in the CS002 local plane
    uvdist = sqrt(sum(uvw_rot.^2, 2) - (uvw_rot * normal).^2); % In meters
	weight = ones (length (uvdist), 1);
	uvec = u(:); vvec = v(:);

	% Find out type of weighting requested
	switch lower (wparm.type)

		case {'uniform'}
			% For every visibility, figure out how many other sampled 
			% visibilities fall within the cell of radius wparm.cellrad
			for vis = 1:length (uvdist)
				% Checking all visibilities, bit bruteforce
				dist = hypot ((uvec(vis)-uvec), (vvec(vis)-vvec));
				weight (vis) = length (find (dist < wparm.cellrad));
			end;
		
		otherwise
			fprintf (2, 'Weighting scheme %s is not known!\n', wparm.type);
	end;

	if (~plt)
		imagesc (reshape (weight, [288 288])); colorbar;		
		title (sprintf ('%s weighting, cellrad: %d', wparm.type, wparm.cellrad));

		plot (uvdist, weight);
		xlabel ('UVdistance (m)'); ylabel ('Weight (no. of visibilities');
		title ('Weights against UV distance');
	end;
