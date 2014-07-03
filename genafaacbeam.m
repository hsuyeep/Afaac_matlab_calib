% Script to generate the average primary beam from AARTFAAC LBAs. The basic input
% is a full EM simulation which incorporates the full aartfaac array geometry,
% both polarizations oriented correctly wrt. each other, an accurate ground plane
% geometry, an approximate soil impedance and mutual coupling between the 
% elements. This is in contrast to the earlier simulation (available in 
% LBA_beam/CS1/calculateLBAbeam.m) which assumed an infinite ground plane, and 
% the field configuration of CS1. The simulation has been carried out in WIPL-D
% by Michel Arts (arts@astron.nl).
%
%  Note that since the simulation was generated at a fixed azi/el + freq grid, 
% interpolated values are returned for the user defined azi/el and freq. 

% Arguments: 
%   l,m : The direction cosine grid on which the primary beam is to be evaluated.
%  freq : The vector of frequencies (Hz) at which the beam will be generated.
% array : The AARTFAAC configuration for which the beam needs to be generated.
%         One of 'lba_outer', 'lba_inner';
% pol   : 'X' to generate an X pol beam, 'Y' for Y pol.
% deb   : Generates debug plots. 1 => final patterns, 2 => per dipole patterns.

% Returns:
% theta  : The theta axis of E_avg11
% phi    : The Phi axis of E_avg11
% E_avg11: The output power at the simulated azi/el + freq, based on hardcoded 
%		   antenna impedance.
%  beam : An lxm matrix containing the power pattern averaged over the chosen 
%         array configuration.
% pep/04Mar14
function [theta, phi, E_avg11, beam] = genafaacbeam (l, m, freq, array, pol, deb)
	persistent simload;
	persistent Positions;
	persistent Freq;
	persistent Phi;
	persistent Theta;
	persistent Zant;
	persistent g_z;

	% Generate the field pattern at values available in the simulation.
	if (isempty(simload))
		fprintf (2, 'Loading simulation from  LBA_core.mat\n');
		load ('LBA_beam/arts_sim/LBA_core.mat'); % Load the simulation products
		simload = 1;
	end;

	% Generate list of dipoles in specified array config.
	switch array
		case 'lba_inner'    % The first 48 of 96 dipoles per station
			sel = [1:2:576];

		case 'lba_outer'
			sel = [577:2:1152];
			
	end;
	
	% Choose polarization. Default above is X-pol
	% The even numbered elements are in Phi=135deg plane, 
	% called the Y pol. here.
	if pol == 'Y'	
		sel = sel + 1;
	end;

	% Generate output voltages from LNAs.
	% Nelem = length (sel); % Number of selected antenna elements.
	Nelem = 1152;
	Glna = 1; % Define Voltage gains of the LNAs
	Vout=zeros(Nelem, 2, length(Theta), length(Phi), length(Freq));
	E11 = zeros (Nelem, length(Theta), length(Phi), length(Freq));
	for NF=1:length(Freq)
		Zlna = (700./(1+i*2*pi*Freq(NF)*700*15e-12))*eye (Nelem); % Define impedances of the LNAs.
		% B=Glna*Zlna/(Zant(sel,sel,NF)+Zlna); % Is it correct to ignore the other cross terms?
		B=Glna*Zlna/(Zant(:,:,NF)+Zlna); % Is it correct to ignore the other cross terms?
		for NT=1:length(Theta)
			for NP=1:length(Phi)
	       		Vout(:,:,NT,NP,NF)=B*g_z(:,:,NT,NP,NF);
				E11 (:,NT,NP,NF) = ... 
				sqrt(abs(Vout (:, 1, NT, NP, NF)).^2 + abs(Vout(:,2,NT,NP,NF)).^2);
			end
		end
		% Should be averaged over only chosen set of antennas.
		E_avg11 = squeeze (mean (E11(sel,:,:,:), 1)); 
	end

	if (deb > 1)
		% Plot individual dipole patterns
		filename = 'dipbeams.gif';
		fprintf ('DEBUG: Generating GIF of individual dipole patterns.\n');
		for ind = 1:size (E11, 1)
			figure(1)
			% subplot (121);
			plot (Positions(1,sel), Positions(2,sel), '.');
			hold on;
			plot (Positions(1,sel(ind)), Positions(2,sel(ind)), 'r.');
			xlabel ('m'); ylabel ('m'); 
		 	title (sprintf ('AARTFAAC %s config.\n', array));
			drawnow;

			% subplot (122);
			figure (2);
		    imagesc (Phi, Theta, 10*log10(squeeze(abs(E11(ind, :,:, 8))))); colorbar;
			xlabel ('Phi (deg)'); ylabel ('Theta(deg)');
			title (sprintf ('Dipole %d of config %s, at Freq %f.', ...
			  			  ind, array, Freq(8)));
		    drawnow;

%		    frame = getframe(1);
%		    im = frame2im(frame);
%		    [imind,cm] = rgb2ind(im,256);
%		    if ind == 1;
%		        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%		    else
%		        imwrite(imind,cm,filename,'gif','WriteMode','append');
%		    end
		end

	end;

	% Generate the azi/el values desired by the user, as specified in l,m
	[lgrid, mgrid] = meshgrid(l, m);
	dist = sqrt(lgrid.^2 + mgrid.^2);
	thetai = asin(dist);
	thetai(dist >= 1) = 0;
	thetai = thetai * 180/pi; % Convert to degrees.
	phii = mod(atan2(lgrid, mgrid), 2 * pi)*180/pi; % Convert to degrees.

	% Interpolate the simulated beam to the frequency and l,m positions needed 
	% by the user.
	beam = zeros(length(l), length(m), length(freq));
	for idx = 1:length(freq)
	    beam(:, :, idx) = interp3 (Phi, Theta, Freq, E_avg11, phii, thetai, freq(idx) * ones(size(phii))) .* (dist < 1);
	end
	
	phi = Phi;
	theta = Theta;

	% Write out generated pattern over frequency as a gif
	if (deb > 0)
		filename = 'powerpatt.gif';
		for ind = 1:length (freq)
			figure(1)
		    imagesc (l, m, 10*log10(beam(:,:,ind)./max(max(beam(:,:,ind)))));colorbar;
			xlabel ('l'); ylabel ('m');
			title (sprintf ('Normalized %s power pattern(dB): %s at %d Hz.', pol, array, freq(ind)));
		    drawnow;

		    frame = getframe(1);
		    im = frame2im(frame);
		    [imind,cm] = rgb2ind(im,256);
		    if ind == 1;
		        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
		    else
		        imwrite(imind,cm,filename,'gif','WriteMode','append');
		    end
		end
	end;
