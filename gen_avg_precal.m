% Function to generate an averaged calibration solution, after carrying out
% Time domain sigma clipping over the window of temporal extent of the complex
% gains.

% Arguments:
% 	gainwindow:  Complex array containing the instantaneous estimated gains 
%				 over a temporal window.
% Returns:
%	avg_cal   :  The per antenna average complex gain over the temporal window.
%	std_cal	  :  The per antenna complex gain st. dev. over the temporal window.
% pep/06Oct12

function [avg_cal, std_cal] = gen_avg_precal (gainwindow)
	ntime = size (gainwindow, 1);
	nant   = size (gainwindow, 2); % cols are antennas, rows are times.
	thresh = 3; % sigma
	maxiter= 5;
	prev_avg = 0; prev_std = 0;
	for i = [1:nant]
		selgains = abs(gainwindow (:, i)); % NOTE: Currently working on amps!
		for j = [1:maxiter]         % sel == 1 => take the value
			
			avg_cal = mean (selgains); % avg. cmplx gain over window for ant i.
			std_cal = std  (selgains); % std. over window for ant i.
			selgains = selgains (selgains < thresh*std_cal + avg_cal);
			disp (['Iter ' num2str(j) ',' num2str(length(selgains)) ...
				   ' gains remaining.']);
			disp (['Avg: ' num2str(avg_cal) 'Std: ' num2str(std_cal)]);
			if (prev_avg == avg_cal) 
				disp (['Breaking after ' num2str(j) ' loops. Avg: ' ... 
						num2str(avg_cal) ' Std:' num2str(std_cal)]);
				break;
			end;
			prev_avg = avg_cal; 
 			prev_std = std_cal;
		end;
	end;
