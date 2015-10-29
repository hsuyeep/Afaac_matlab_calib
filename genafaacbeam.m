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

% Details on the simulation and how to obtain the beam pattern from it.
% ----------------------------------------------------------------
% - Simulation output is an open circuit voltage, i.e., voltage induced in the
% antenna by a point source with a field of 1V/m located at each of the probed
% positions.
% - To calculate the beam pattern, we need to calculate the output voltage of
% each dipole after the LNA. This is done by incorporating the gain of each LNA (set to 1) and impedance
% of each antenna (available in the Zant matrix) into a combined gain and multiplying the opencircuit voltage with it.
% This results in a gain matrix G, i.e., G = B.g_z;
% 
% - The basic equation:
%     V = G.E, 
%     where V = [v_x v_y] are the two polarizations of the received voltage
%     from the output of the antenna.
%     G = [ g^theta_x g^phi_x
%           g^theta_y g^phi_y] are the gains per theta or phi component and X or Y polarization.            
%     E = [E_theta E_phi] are the components of the incoming electric field,
%          resolved in directions perpendicular to directions of propagation,
% 
% - To generate the Stokes components, we need the Coherency Matrix, defined as 
%   V.V^H = G.E.E^H.G^H = [v_x.v_x* v_x.v_y*
%                          v_y.v_x* v_y.v_y*];
%   For a stokes-I beam pattern,  E.E^H can be replaced by
%                                           [0.5  0;
%                                            0  0.5], as then the we get
%                            stokes-I = v_x.v_x* + v_y.v_y*
%
%  
% --------------------------------------------------------------------
% Michel Arts's email of 21Oct2013 wrt. AARTFAAC-6 beam response simulation:
% Hello Peeyush and Stefan,
%  
% I have put the simulation data on the ASTRON FTP-server. You can
% download the mat-file via the link
%  
% ftp://ftp.astron.nl/outgoing/arts/LBA_core.mat
%  
% The file LBA_core.mat contains the following variables:
%  
% Freq: vector with the simulated frequencies in Hz.
%  
% Phi: the simulated azimuth angles in degrees. With respect to the array
% configuration the positive x-axis points to the east and the positive
% y-axis points to the north. The angle Phi is defined w.r.t. the x-axis
% as in a standard spherical coordinate system.
%  
% Theta: the simulated elevation angles in degrees (w.r.t. zenith).
%  
% Positions: coordinates of the antenna elements w.r.t. element 1 in
% meter. There are 1152 elements (576 LBA's each containing 2 elements).
% The odd numbered elements are in the plane Phi=45 degrees. The even
% numbered elements are in the plane Phi=135 degrees. The angle Phi is the
% angle w.r.t. the x-axis as defined in a standard spherical coordinate
% system. Element n and n+1 with n an odd number, have the same
% coordinates because these elements are at the same position. The only
% difference is that antenna n+1 is rotated 90 degrees w.r.t. antenna n.
%  
% Zant: the impedance matrix of the antenna. The format is element index
% x element index x frequency index.
%  
% g_z: the open circuit voltages of each element for an incoming plane
% wave with an electric field strength of 1 V/m. The format is element
% index x component index x theta index x phi index x frequency index.
% Component index=1 corresponds with the theta-component, component
% index=2 corresponds with the phi-component.
%  
% You can calculate the output voltage of each LNA by the following
% code:
%  
% Vout=zeros(size(Zant,1),2,length(Theta),length(Phi),length(Freq));
% for NF=1:length(Freq)
%    B=Glna*Zlna/(Zant(:,:,NF)+Zlna);
%    for NT=1:length(Theta)
%        for NP=1:length(Phi)
%            Vout(:,:,NT,NP,NF)=B*g_z(:,:,NT,NP,NF);
%        end
%    end
% end
%  
% Zlna is a diagonal matrix with the LNA-impedances on the diagonal. Glna
% is the voltage gain of the LNA (scalar). The auxiliary variable B is
% used to speed up the calculation by minimizing the number of
% multiplications in the inner for-loop.  Vout(:,1,NT,NP,NF) contains the
% output voltages of the LNAs if the array is excited by a planar wave at
% frequency Freq(NF) incident at angle Theta(NT),Phi(NP) with a
% theta-componant of 1 V/m and a phi-component of 0 V/m..
% Vout(:,2,NT,NP,NF) contains the output voltages of the LNAs if the array
% is excited by a planar wave at frequency Freq(NF) incident at angle
% Theta(NT),Phi(NP) with a theta-componant of 0 V/m and a phi-component of
% 1 V/m..
%  
% Michel
% -------------------------------------------------------------------------
% Followup on sharp change in beamshape from 50 to 55 MHz. Resulted in a finer
% simulation  between 40 and 60 MHz.
% Hi Peeyush,
%  
%  Just wondering if you have any thoughts on the sudden change in the
% beam
% properties at 55 MHz in the new simulations?
%  
%  
% A possibility is the finite ground plane. It is 3 x 3 meter. This is
% half a wavelength at 50 MHz. So it is likely that the change in
% radiation pattern is caused by resonance effect of the ground plane.
%  
% Michel

function [theta, phi, E_avg11, beam, v_ant_lna_avg] = genafaacbeam (l, m, freq, array, pol, deb)
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

	% Generate output voltages from LNAs.
	% Nelem = length (sel); % Number of selected antenna elements.
	Nelem = 1152;
	Glna = 1; % Define Voltage gains of the LNAs
	Vout=zeros(Nelem, 2, length(Theta), length(Phi), length(Freq));
    G = zeros (length(Theta), length(Phi), length(Freq), Nelem/2, 2, 2);
    
    % Generate Coherency matrix and stokes-I beam for all antennas
    CM = zeros (size (G));
    stI = zeros (Nelem/2, length (Theta), length(Phi), length(Freq)); 
    
	% E11 = zeros (Nelem, length(Theta), length(Phi), length(Freq));
	for NF=1:length(Freq)
        % Zlna taken from the Master's thesis of Maria Krause, pg. 27
        % (Schematic) and pg. 47 for estimated values of R and C.
		Zlna = (700./(1+i*2*pi*Freq(NF)*700*15e-12))*eye (Nelem); % Define impedances of the LNAs.
		
		B=Glna*Zlna/(Zant(:,:,NF)+Zlna); % Is it correct to ignore the other cross terms?
		for NT=1:length(Theta)
			for NP=1:length(Phi)
	       		Vout(:,:,NT,NP,NF)=B*g_z(:,:,NT,NP,NF);
				
			end
        endfprintf (1, '<-- Generated output voltages per dipole, per E-field component.\n');
    % Generate Gain matrix (see header comment in genafaacbeam.m) per
    % antenna, for all antennas and both pols.
    
    
    for NF=1:length(Freq) 
        for NT=1:length(Theta)
            for NP=1:length(Phi)
                for ant=1:2:Nelem % Because we use up both pols of the ant.
                G(NT, NP, NF, (ant+1)/2, :, :) = 
                 [squeeze(Vout(ant  ,1,NT,NP,NF)) squeeze(Vout(ant  ,2,NT,NP,NF));
                  squeeze(Vout(ant+1,1,NT,NP,NF)) squeeze(Vout(ant+1,2,NT,NP,NF))];
                end;
            end;
        end;
    end;
    fprintf (1, '<-- Calculated full pol. Jones gain matrix per antenna.\n');
    
       
    for NF = 1:length(Freq)
        for NT=1:length(Theta)
            for NP=1:length(Phi)
                for ant = 1:Nelem/2
                    tmp = squeeze(G(NT,NP,NF,ant,:,:)); % 2x2 Jones matrix
                    CM(NT,NP,NF,ant,:,:) = tmp*0.5*eye(2)*tmp';
                    % below is already a power due to adding up xx and yy
                    % of the CM.
                    stI (ant,NT, NP, NF) = CM(NT,NP,NF,ant,1,1) + CM(NT,NP,NF,ant,2,2); 
                end;
            end;
        end;
    end;
    fprintf (1, '<-- Generated coherency matrix and stokes-I power beam.\n');
        
      
      
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
