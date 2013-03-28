% Script to generate a 2D gaussian with the specified parameters.
% pep/20Dec12

%    """Return a 2D Gaussian function with the given parameters.
%
%    Args:
%
%        height (float): (z-)value of the 2D Gaussian
%
%        center_x (float): x center of the Gaussian
%
%        center_y (float): y center of the Gaussian
%
%        semimajor (float): major axis of the Gaussian
%
%        semiminor (float): minor axis of the Gaussian
%
%        theta (float): angle of the 2D Gaussian in radians, measured
%            between the semi-major and y axes, in counterclockwise
%            direction.
%
%    Returns:
%
%        (function): 2D Gaussian
%    """

function [gau] = gaussian(x, y, height, center_x, center_y, semimajor, semiminor, theta)
	% NOTE: Put a check on the major/minor axis actually exceeding the extent
	% of X and Y provided, and scream in that case.
	gau = zeros (length(x), length(y));
	
	for iind = 1:length(x)
		for jind = 1:length(y)
			gau(iind, jind) = height * exp( -log(2.0) * (((cos(theta) * (x(iind) - center_x) + sin(theta) * (y(jind) - center_y)) / semiminor)^2.0 + ((cos(theta) * (y(jind) - center_y) - sin(theta) * (x(iind) - center_x)) / semimajor)^2.));
		end;
	end;
