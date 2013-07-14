% Script to fit a 2D gaussian on a portion of the image containing
% sources. Basic fit: no position angle parameter estimated.
% pep/21Dec12

% Arguments:
%	image:	2D array containing source islands (noise?)
%	    l:  l-coordinates of this section of the image
%	    m:  m-coordinates of this section of the image
% 	init_par: The initial estimates of the parameters. 
%				Suitably initialized.
% init_par = [peak_x, peak_y, peak_flux, sigx, sigy], where:
%	peak_x = x coordinate of peak pixel, used for putting offsets to model.
%	peak_y = y coordinate of peak pixel, used for putting offsets to model.
%	peak_flux = flux of peak pixel in subimage.
%	sigx = spread along x-axis (across columns of an image).
%	sigy = spread along y-axis (across rows of an image).


% Returns:
%	fitparms: A structure containing an estimate of the same parameters as 
%			init_par, returned in the same order. 
%     resid : Frobenius norm of the residual image, for fit quality check.

% NOTE: Can improve fit parameters by including a position angle in the model;
% see gaussian.m. Also, by providing initial estimates via moments.
% NOTE: How do we know if the fit did not go well?

function [fitparams, resid] = fit2dgauss(image, l, m, init_par, debug,debwinhdl)

% Debug function, uncomment all commented sections, then run as fit2dgauss()
%function [fitparams] = fit2dgauss() 
%	l = -1:.2:1;
%	m = -1:.2:1;
%	debug = 1;

	% NOTE: Assumes l, m are vectors.
	[lmesh, mmesh] = meshgrid (l, m);

%	% Define a specific 2D function with 3 parameters to test "fminsearch"
%	A = 5;
%	sigx = 0.8;
%	sigy = 0.4;
%	x0 = 0; y0 = 0;
%	noise = A/10*randn(size (lmesh)); % Gaussian noise.
%	image = A*exp(-(lmesh - x0).^2/(2*sigx^2)).*exp(-(mmesh-y0).^2/(2*sigy^2)) + noise;
%	
%	init_par = [1 1 1 1 1];
%	imagesc (l, m, image); 
%	figure; 
%	subplot (1,2,1); plot (image (:,int32(size(image,2)/2)+x0));
%	subplot (1,2,2); plot (image (int32(size(image,2)/2)+y0,:));
	

	% Attempt to invoke fminsearch
	fitparams = fminsearch(@residual,init_par, [], lmesh, mmesh, image);
	
	 % For debug:
	model_image = init_par(3) * ...
					exp(-(lmesh-init_par(1)).^2/(2*init_par(4)^2)).*...
					exp(-(mmesh-init_par(2)).^2/(2*init_par(5)^2));
	fit_image = fitparams(3) * ...
					exp(-(lmesh-fitparams(1)).^2/(2*fitparams(4)^2)).*...
					exp(-(mmesh-fitparams(2)).^2/(2*fitparams(5)^2));
	resid_image = image - fit_image;
	resid = norm (resid_image, 'fro');
	if debug > 1
		figure (debwinhdl); 
		subplot (131); imagesc (image); title ('Data'); colorbar;
		subplot (132); imagesc (fit_image); title ('Fit image'); colorbar;
		subplot (133); 
		imagesc (resid_image); title ('Residual image'); colorbar;
	 end;

	% Define a residual cost function for minimization.
	function err = residual(par, l, m, image)
	  % global X Y Zdata;
		model_image = par(3) * ...
						exp(-(l-par(1)).^2/(2*par(4)^2)) .* ...
						exp(-(m-par(2)).^2/(2*par(5)^2));
		err = sum(sum( (model_image - image).^2 ));
