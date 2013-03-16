% Script to flag dead visibilities, based on average levels over all 
% visibilities in a single timeslice.
% Arguments:
%	acc : ACM for this time instant
%	tobs: Time of obs
% 	freq: Freq. of obs.
% thresh: Threshold of median cutoff. All vis>thresh*median are flagged.

% Returns:
%	uvflag: flag matrix with 1s indicating visibilities to flag.
%	flagant: indices of antennas reported to be bad.
% pep/30Jan13

function [uvflag, flagant] = flagdeadcorr (acc, tobs, freq, hithresh, lothresh)

	ac = abs (acc - diag(diag(acc)));
	uvflag = zeros (size (acc)); % 1 => bad visibility
	flagant = [1:length(acc)];

	% Find missing antennas. Row 20 is an arbit choice.
%	missant = ac (20, :) == 0; 
%	missant(20) = 0; % Since missing autocorr will also be marked as bad.

	% Find antennas with very low values.
	mu = mean (ac); med = median (mu);
	fl = mu < lothresh*med;
	% missant = missant | fl;

	flagant = flagant(fl);

 
	uvflag (fl == 1, :) = 1; % Update uvflag with missing ants.
	uvflag (:, fl == 1) = 1;

	% For each antenna, find median absolute values of visibilities with all
	% other antennas.
	for ind = 1:size (acc, 1)
		ant = ac (ind, :);
		medant = median (ant); % NOTE: Median is unbiased wrt. outliers.
		% plot (ant);
		uvflag (ind, ant > medant + hithresh*medant) = 1; % Flag discrepant vis.
		uvflag (ind, ant < medant - hithresh*medant) = 1;
	end;
