% Function to obtain the noise floor and dynamic range of an image. The stats.
% are generated from a central region in the image (with a specifyable size), 
% and by carrying out sigma clipping (specifyable sigma) to root out possible 
% sources which can bias the statistics.
% Arguments:
%	map     : mxm map, usually the output of an imaging task
% 	cellsize: halfside of a square region for estimating the variance within the 
%			  image.
%   thresh  : Sigma clipping threshold.
% pep/24Oct12

function [dr, sig] = getimagedr (map, cellsize, thresh)
	% Search within the center of the image, as noise is expected to be better 
	% behaved there.	
	maxpix = max (map(:));
	cen = [size(map)/2];
	tmp = map (cen(1)-cellsize:cen(1)+cellsize, cen(2)-cellsize:cen(2)+cellsize);
	reg = tmp(:);
	mask = ones (size (reg));
	
	for ind = 1:5
		sig = std (reg(mask == 1));
		avg = mean (reg (mask == 1));
		mask = mask & (reg < avg + thresh*sig);
	end
	dr = maxpix/sig;
