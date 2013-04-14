% Script to establish the goodness of a calibration solution, based on 
% comparisons against a sliding window of good solutions.
% pep/25Mar13
% Arguments:
%   solwindow: The solution window consisting of an array of solution records
%   currsol  : The current calibration solution on which judgement needs to be 
%				passed.
%  gainmask  : Mask holding flagged antennas. Not passed through solparm due to 
%			   duplication of gainmask several times. 
%   solparm  : Structure containing the thresholds used to determine goodness.
%
%  solparm.solthresh: Threshold on the std. of the variation of gain solutions 
%			   across a station. Assumed that a station is mostly calibrated.
%  solparm.errthresh: Threshold on the relative error on current calsol. wrt. 
% 			   mean calsol.
%   debug  : Debug level of script. Uses handles passed in solparm for plotting 
%			 thresholds against data.
% Returns:
% solstat.goodstncal:Bool desc. the goodness of a calsol wrt.intra station test.
% solstat.err       : Relative error in percent. wrt. mean solution over window.
% solstat.stat_std_arr: Complex stds. of each station over window.
% solstat.mse_amp   : MSE over amplitudes
% solstat.mse_ph    : MSE over phases
% solstat.mse_re    : MSE over real components
% solstat.mse_im    : MSE over imag components
% solstat.mse       : MSE over cartesian differences.
% Algorithm:
%   The check for a good solution occurs on the phase and amplitude of the 
% solution. First, the relative error between the current solution and the mean
% solution is determined. If this is large, the phase is checked with a more 
% stringent condition (since a phase error is more damaging than an amp. error).

function [solstat] =  goodcalsol (solwindow, nsol, currsol, gainmask, solparm, debug)
	
	if (nsol == 0) % Window is empty
		% NOTE: Order in which structure entries are filled matters! CRAZY!
		solstat.mse = 0;
		% Currently not being u
		solstat.mse_re  = 0;
		solstat.mse_im  = 0;
		solstat.mse_amp = 0; 
		solstat.mse_ph  = 0; 
		solstat.err = 0;
		solstat.pherr = 0
		solstat.stat_std_arr = zeros (1, 6);
		% No window to compare against, just check the consistency over stations
		% of variance in gain amp and phase within a station.
		for sta=0:5
			sta_ind = [sta*48+1 : (sta+1)*48];
			% Use only unflagged antennas per station, check phase variance.
			% Take care of phase wraps.
			station_ph = angle (currsol (gainmask(sta_ind) == 0));
			station_amp = abs (currsol (gainmask (sta_ind) == 0));
			
			% Phase unwrap
			diffphfact = abs(station_ph ./ pi);
			sel = diffphfact > 0.5; % Choose phase differences > 0.5 rad.
			diffphfact (sel) = ceil (diffphfact (sel)); % Convert to integers.
			diffphfact (sel == 0) = floor (diffphfact(sel == 0));
			% Choose +ve phase differences, need to subtract 2pi mult.
			sel = station_ph > 0; 
			station_ph (sel) = station_ph (sel) - diffphfact(sel)*pi;
			station_ph (sel==0) = station_ph (sel==0) + diffphfact(sel==0)*pi;
			station_ph = station_ph * (180/pi); % Convert to degrees.

			solstat.stat_std_arr(sta+1) = std(station_ph);
	
		end;

		% Why should the std be comparable to the mse threshold?
		% NOte this can be used to eliminate the sols. with completely random 
		% phases, especially when no window is available...
		if (sum (solstat.stat_std_arr > solparm.mse_ph_thresh)>	1)
			solstat.goodstncal = 0;
		else 
			solstat.goodstncal = 1;
		end;
		return;
	end;


	% Compute mean of real and imaginary components separately.
	if (nsol == 1)
		meansol = solwindow(1,:);
	else 
		meansol = mean(solwindow (1:nsol,:));
	end;

%	if (debug >= 4)
%		figure (solparm.gainplt);
%		subplot (1, 2, 1);
%		plot (abs(solwindow)');
%		subplot (1, 2, 2);	
%		plot (angle (solwindow)');
%
%		figure (solparm.currsolplt);
%		subplot (1,2,1);
%		plot (abs(meansol), '-b');
%		hold on;
%		plot (abs(currsol), '-r');
%		hold off;
%		subplot (1,2,2);
%		plot (angle(meansol), '-b');
%		hold on;
%		plot (angle(currsol), '-r');
%		hold off;
%	end;

	% Generate mean squared error (mse) and sum of squared error (sse).
	diffcart = meansol - currsol; % Cartesian complex difference 
	diffamp = abs (meansol) - abs (currsol); % Ampl. difference.
	diffph  = (angle (meansol) - angle (currsol)); 

	% Take care of phase wraps.
	diffphfact = abs(diffph ./ pi);
	sel = diffphfact > 0.5; % Choose phase differences > 0.5 rad.
	diffphfact (sel) = ceil (diffphfact (sel)); % Convert to integers.
	diffphfact (sel == 0) = floor (diffphfact(sel == 0));
	sel = diffph > 0; % Choose +ve phase differences, need to subtract 2pi mult.
	diffph (sel) = diffph (sel) - diffphfact(sel)*pi;
	diffph (sel == 0) = diffph (sel == 0) + diffphfact(sel == 0)*pi;
	diffph = diffph * (180/pi); % Convert to degrees.
	[m, v, sel] = robustmean (diffph, 5); % Ignoring phases > 5sigma.

	% Arbit. declaring bad if >0.5 of antennas need to be thrown while 
	% constructing a robust mean.
%	if (sum(sel) < length (diffph)/2)     
%	end;

	npar = length(meansol); 

	% 2 for re/im, avg. err per ant.
	solstat.mse = (diffcart * diffcart')/(2*npar); 
	% Currently not being used.   
	solstat.mse_re  = sum (real (diffcart).^2)/npar; % avg. err in real 
	solstat.mse_im  = sum (imag (diffcart).^2)/npar; % avg. err in imag 

	solstat.mse_amp = sum (diffamp.^2)/npar;         % avg. err in amp.
	solstat.mse_ph  = sum (diffph(sel == 1).^2)/npar;% avg. err in ph, in deg.

	% Compute relative error of this solution wrt. mean solution.
	solstat.err = 100*sum(abs(meansol - currsol)) / sum(abs(meansol));
	% If phases are constant, but amplitudes change, it is still a good sol.
	solstat.pherr = 100*sum(angle(meansol - currsol)) / sum(abs(meansol));

	% Compute per station gain variances. NOTE: All solution always have 288 
	% elements. Also, std (complex number) = sqrt (std(re).^2 + std(im).^2);
	for sta=0:5
		sta_ind = [sta*48+1 : (sta+1)*48];
		% Use only unflagged antennas per station
		stat_std = std(currsol (gainmask(sta_ind) == 0));
		stat_std_re = std(real(currsol (gainmask(sta_ind) == 0)));
		stat_std_im = std(imag(currsol (gainmask(sta_ind) == 0)));

		solwin_std = std (meansol (gainmask(sta_ind) == 0));
		solwin_std_re = std(real (meansol(gainmask(sta_ind) == 0)));
		solwin_std_im = std(imag (meansol(gainmask(sta_ind) == 0)));

		solstat.stat_std_arr(sta+1) = stat_std;

		% Print the offending stds.
%		fprintf (1, '%3.1f/%3.1f [%3.1f/%3.1f, %3.1f/%3.1f] ', ... 
%				stat_std, solwin_std, stat_std_re, solwin_std_re, ...
%				stat_std_re, solwin_std_re, stat_std_im, solwin_std_im);

		% NOTE: comparison is true only if all vector components are true.
		% if ((stat_std  < solwin_std-solparm.solthresh*solwin_std) ||  ...
		%     (stat_std  > solwin_std+solparm.solthresh*solwin_std)) 
%		if (abs(1 - stat_std/solwin_std) > solparm.solthresh)
%			goodstncal = 0;
%			% pause;
%
%			%stat_std_arr (:)= 0; % stn std should be ignored if goodstncal == 0
%			% break;
%		end;
		% solwin_std_arr(sta+1) = solwin_std;
	end;
	
	% Also considered a bad solution if relative error is larger than specified.
	solstat.goodstncal = 1;
	if (solstat.goodstncal == 1 && solstat.err > solparm.errthresh)
		if (solstat.mse_ph < solparm.mse_ph_thresh) % < 30deg. phase err is OK to pass.
			solstat.goodstncal = 1;
		else
			solstat.goodstncal = 0;
		end;
		% pause;
	end;
	fprintf (1, '{e:%5.1f, %4.1f, m:%5.2f, r:%5.2f, i:%5.2f, a:%5.2f, p:%5.2f} ', solstat.err, solstat.pherr, solstat.mse, solstat.mse_re, solstat.mse_im, solstat.mse_amp, solstat.mse_ph);

	% Uncomment conditional to plot only the bad solutions.
	if (debug >= 4) % && goodstncal == 0)
		
		figure (solparm.gainplt);
		subplot (1, 2, 1);
		% plot (abs(solwindow)');
		plot (real(solwindow)');
		xlabel ('Antenna number');
		title ('Real component of sol window');
		subplot (1, 2, 2);	
		% plot (angle (solwindow)');
		plot (imag (solwindow)');
		xlabel ('Antenna number');
		title ('Imag component of sol window');
		
		figure (solparm.currsolplt);
		subplot (1,2,1);
		plot (abs(meansol), '-b');
		xlabel ('Antenna number');
		title ('Mean solution Vs. current solution');
		hold on;
		plot (abs(currsol), '-r');
		hold off;
		subplot (1,2,2);
		plot (angle(meansol), '-b');
		hold on;
		plot (angle(currsol), '-r');
		hold off;
		drawnow;
		% pause;
	end;
