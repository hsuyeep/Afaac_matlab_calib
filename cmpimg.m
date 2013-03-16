% Script to generate relative image statistics from provided calibrated ACMs.
% Note that calibration quality is better compared using calibration solutions,
% rather than images, which average over all gain solutions to generate every
% output pixel value.
% pep/08Mar13
%
% Arguments:
%  acc  : multidimensional matrix of ACMs to compare
%  nacc : The number of provided accs.
%imgparm: Structure containing imaging parameters. If passed empty,
%  		  default values are assigned.
%  tobs : Time of observation in mjdsec
%  freq : Frequency of observation in Hz.
%flagant: list of flagged antennas, assumed identical for all acc. Used only if
% 		  imgparm is not provided.
%
% Returns:
%  img.DR     : An array containing the dynamic range of all images wrt. 1st
%  img.diffmap: Multidimensional matrix returning the difference of all images wrt.
%		  1st.
% 
function [img, imgparm] = cmpimg (acc, nacc, imgparm, tobs, freq, flagant)

	lambda = 299792458./freq; 		% in m.

	% Check if imgparm is empty, if yes, generate default.
	if (isempty(imgparm) == 1)
		imgparm.dofft = 1;    % Carry out FFT imaging
    	imgparm.duv=lambda./2;% Default image just the full Fov (-1<l<1)
    	imgparm.Nuv = 500;    % size of gridded visibility matrix
    	imgparm.uvpad = 512;  % specifies if any padding needs to be added

		% dimensionless, in dir. cos. units
		if (imgparm.dofft == 1)
	    	dl = 299792458./(freq .* (imgparm.uvpad .* imgparm.duv)); 
			img.map = zeros (imgparm.uvpad, imgparm.uvpad, nacc);
			img.diffmap = img.map;
		else	
			img.l = linspace (-1, 1, 512);
			img.m = img.l;
			img.map = zeros (length(img.l), length(img.m), nacc);
			img.diffmap = img.map;
		end;
		
		% Assuming LBA_OUTER array
    	load ('poslocal.mat', 'posITRF', 'poslocal'); 
		if (imgparm.dofft == 1)
			uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
			vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
			% Generate flagged positions
			[imgparm.uloc_fl, imgparm.vloc_fl] = ... 
					gen_flagged_uvloc (uloc(:), vloc(:), flagant); 
		else
			sel = ones (288, 1);
			sel (flagant) = 0;
			imgparm.uloc_fl = poslocal (sel>0,1);
			imgparm.vloc_fl = poslocal (sel>0,2);
		end;
		imgparm.radec = 0; % Only zenith maps
	end;

	if (imgparm.dofft == 1)
		fprintf (2, 'Forming FFT images!\n');
	else
		fprintf (2, 'Forming DFT images!\n');
	end;

	% Form images from all the ACMs.
	for ind = 1:nacc
		if (imgparm.dofft == 1)
   			[radecmap, img.map(:,:,ind), calvis, img.l, img.m] = ... 
		 fft_imager_sjw_radec (acc(:,:,ind),imgparm.uloc_fl(:), ... 
			imgparm.vloc_fl(:), imgparm.duv(ind), imgparm.Nuv, imgparm.uvpad,...
				tobs(ind), freq(ind), imgparm.radec);
			img.map(isnan(img.map) == 1) = 0;
		else
			img.map(:,:,ind) = acm2skyimage (acc(:,:,ind),imgparm.uloc_fl, ... 
				imgparm.vloc_fl, freq(ind), img.l, img.m);	
			img.map(isnan(img.map) == 1) = 0;
		end;
	end;	
 
	% Generate statistics on images
	for ind = 2:nacc
		% Diffmap
		img.diffmap (:,:, ind) = abs(img.map(:,:,ind) - img.map(:,:,ind-1)).^2;

		% Normalized Mean square error estimate = normalised least squares err
		img.nlse(ind) = ... 
			sum (sum (abs(img.map(:,:,ind) - img.map(:,:,ind-1)).^2)) / ...
			sum (sum (abs(img.map(:,:,ind-1)).^2));
							 
	end;	
