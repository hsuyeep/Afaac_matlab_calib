% Script to generate the beam for the AARTFAAC zenith pointing.
% Beam generated in (l,m) coordinates, based on different antenna 
% configurations.
% pep/20Dec12.

% Arguments:
%	arrayconfig: one of 'LBA_OUTER', 'LBA_INNER', ... 
%	freq	   : Frequency of observation (Hz)
%	grid	   : 'None', 
%	Nuv		   : Nuv (Number of grid points after grid convoution)

% Returns:
%	beam	   : 2D matrix with psf, with peak normalized to 1, and beam power
%				 in linear units.

function [beam] = aartfaac_beam (arrayconfig, grid, freq, Nuv)
	C = 299792458; 
	lambda = C/freq;	% In m.
	duv = 0.5; % lambda/2; 	% In m, Always Nyquist sample the UV plane.

	% Find out array configuration
	if arrayconfig == 'LBA_OUTER'
		fid = fopen ('LBA_OUTER_config.txt', 'rt');
		posITRF = textscan (fid, '%s = [%14.6f64,%14.6f64,%14.6f64]', ...
							'CommentStyle', '#');
	elseif arrayconfig == 'LBA_INNER'
		fid = fopen ('LBA_OUTER_config.txt', 'rt');
		posITRF = textscan (fid, '%s = [%14.6f64,%14.6f64,%14.6f64]', ...
							'CommentStyle', '#');
	end;
	fclose (fid);

	% convert positions to local horizon frame
	% rotation matrix taken from AntennaField.conf file from CS002
	rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
   			   0.9928230000, -0.0954190000, 0.0720990000; ...
           	   0.0000330000,  0.6030780000, 0.7976820000];
	poslocal = [posITRF{2} posITRF{3} posITRF{4}] * rotmat;

	% Generate baseline distances, in lambda units
	uloc = (meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).')/lambda;
	vloc = (meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).')/lambda;
	wloc = (meshgrid (poslocal(:,3)) - meshgrid (poslocal(:,3)).')/lambda;


	% Carry out appropriate gridding and tapering as instructed
	% NOTE: Currently carrying out the weighted nearest neighbour
	% interpolation of Wijnholds.
	vis = zeros (Nuv);
	
    for idx = 1:length(uloc(:))            % For every recorded visibility
    
    	% Get amp. and direction vector of visibility
        ampl = 1; % abs(acc(idx));			
        phasor = 0; % acc(idx) / ampl;
    
    	% Determine the grid along U-axis in which observed visibility falls.
        uidx = uloc(idx) / duv + Nuv / 2;  
        uidxl = floor(uidx);	% Find the lower and higher gridded U-value.
        uidxh = ceil(uidx);
    
        % Find absolute distance of measured visibility from grid points, 
		% in the U-direction only.
        dul = abs(uidx - uidxl);		
        duh = abs(uidx - uidxh);
    
        % Distribute the visiblity amplitude among the two grid points on the 
		% U-axis in proportion to their distance from the observed visiblity.
        sul = duh * ampl;
        suh = dul * ampl;
        
    	% Determine the grid along V-axis in which observed visibility falls.
        vidx = vloc(idx) / duv + Nuv / 2;
        vidxl = floor(vidx);	% Find the lower and higher gridded V-value.
        vidxh = ceil(vidx);
    
        % Find absolute distance of measured visibility from grid points, 
		% in the V-direction only.
        dvl = abs(vidx - vidxl);
        dvh = abs(vidx - vidxh);
    
    	% Distribute the HIGHER u-grid point's share of the observed 
		% visibility amp. between the higher and lower V-grid point.
        sull = dvh * sul;
        suhl = dvh * suh;
    
    	% Distribute the LOWER u-grid point's share of the observed 
		% visibility amp. between the higher and lower V-grid point.
        sulh = dvl * sul;
        suhh = dvl * suh;
        
    	% Now that the observed visiblity amplitude is distributed among its 
    	% surrounding 4 grid points, fill the gridded visibility matrix with
    	% vectors with the same phase as the original observed visibility.
    	% NOTE: Adding the 4 vectors at the corners of the grid square will 
		% give back the original ungridded observed visibility.
    	% NOTE: We need to accumulate this to prevent overwriting the gridded 
		% values from a visibility from a neighbouring grid square.
        vis(uidxl, vidxl) = vis(uidxl, vidxl) + sull; % * phasor;
        vis(uidxl, vidxh) = vis(uidxl, vidxh) + sulh; % * phasor;
        vis(uidxh, vidxl) = vis(uidxh, vidxl) + suhl; % * phasor;
        vis(uidxh, vidxh) = vis(uidxh, vidxh) + suhh; % * phasor;
        
        %W(uidx, vidx) = W(uidx, vidx) + 1;
    end

	% Carry out 2D FFT to generate beam pattern 	
    % FFT imaging: Just take FFT of the gridded visibilities.
	figure;
	subplot (1,2,1);
	imagesc (vis); colorbar;
    vis = conj(flipud(fliplr(fftshift(vis))));
    beam = fftshift(fft2(vis));
	subplot (1,2,2);
	imagesc (20*log10(abs(beam))); colorbar;
