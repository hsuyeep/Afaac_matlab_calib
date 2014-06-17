% Script to generate the average primary beam from AARTFAAC LBAs. The basic input
% is a full EM simulation which incorporates the full aartfaac array geometry,
% both polarizations oriented correctly wrt. each other, an accurate ground plane
% geometry, an approximate soil impedance and mutual coupling between the 
% elements. This is in contrast to the earlier simulation (available in 
% LBA_beam/CS1/calculateLBAbeam.m) which assumed an infinite ground plane, and 
% the field configuration of CS1.
%
%  Note that since the simulation was generated at a fixed azi/el + freq grid, 
% the values returned are after interpolation.
% Arguments: 
%   l,m : The direction cosine grid on which the primary beam is to be evaluated.
%  freq : The Freq (Hz) at which the beam will be generated.
% array : The AARTFAAC configuration for which the beam needs to be generated.
%         One of 'lba_outer', 'lba_inner';
% Returns:
%  beam : An lxm matrix containing the power pattern averaged over the chosen 
%         array configuration.
% pep/04Mar14
function [beam] = genafaacbeam (l, m, freq, array)
	persistent simload;

	% Generate the field pattern at values available in the simulation.
	if (isempty(simload))
		load ('LBA_beam/arts_sim/LBA_core.mat'); % Load the simulation products

		% Generate list of dipoles in specified array config.
		switch array
			case 'lba_inner'
			case 'lba_outer'
				
		end;
	
		Glna = 1; % Define Voltage gains of the LNAs
		Zlna = eye (1152); % Define impedances of the LNAs.
		Vout=zeros(size(Zant,1), 2, length(Theta), length(Phi), length(Freq));
		E11 = zeros (size (Zant,1), length(Theta), length(Phi), length(Freq));
		for NF=1:length(Freq)
			B=Glna*Zlna/(Zant(:,:,NF)+Zlna);
	   		for NT=1:length(Theta)
				for NP=1:length(Phi)
	           		Vout(:,:,NT,NP,NF)=B*g_z(:,:,NT,NP,NF);
					E11 (:,NT,NP,NF) = 0.5 * ...
					(Vout (:, 1, NT, NP, NF) .* conj (Vout(:,1,NT,NP,NF))+...
					 Vout (:, 2, NT, NP, NF) .* conj (Vout(:,2,NT,NP,NF)));
				end
			end
			% Should be averaged over only chosen set of antennas.
			E_avg11 = squeeze (mean (E11, 1)); 
		end
		simload = 1;
	end;

	% Generate the azi/el values desired by the user, as specified in l,m
	[lgrid, mgrid] = meshgrid(l, m);
	dist = sqrt(lgrid.^2 + mgrid.^2);
	thetai = asin(dist);
	thetai(dist >= 1) = 0;
	phii = mod(atan2(lgrid, mgrid), 2 * pi);

	% Interpolate the simulated beam to the frequency and l,m positions needed 
	% by the user.
	beam = zeros(length(l), length(m), length(freq));
	for idx = 1:length(freq)
	    beam(:, :, idx) = interp3 (Phi, Theta, Freq, E_avg11, phii, thetai, freq(idx) * ones(size(phii))) .* (dist < 1);
	end
	
