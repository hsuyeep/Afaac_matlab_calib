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
%     goodstncal : Bool desc. the goodness of a calsol wrt. intra station test.
%     err        : Relative error wrt. mean solution over the window.
%     stat_std   : Complex stds. of each station over window.

function [goodstncal, err, stat_std_arr] =  goodcalsol (solwindow, currsol, ... 
													gainmask, solparm, debug)
	
	% Compute mean of real and imaginary components separately.
	meansol = mean(solwindow);

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

	npar = length(meansol); 
	mse = (diffcart * diffcart')/(2*npar); % 2 for re/im, avg. err per ant.
	% Currently not being used.   
	mse_re  = sum (real (diffcart).^2)/npar; % avg. err in real component
	mse_im  = sum (imag (diffcart).^2)/npar; % avg. err in imag component

	mse_amp = sum (diffamp.^2)/npar;         % avg. err in amp.
	mse_ph  = sum (diffph.^2)/npar;          % avg. err in ph, in deg.

	% Compute relative error of this solution wrt. mean solution.
	err = 100*sum(abs(meansol - currsol)) / sum(abs(meansol));
	% If phases are constant, but amplitudes change, it is still a good sol.
	pherr = 100*sum(angle(meansol - currsol)) / sum(abs(meansol));

	% Compute per station gain variances. NOTE: All solution always have 288 
	% elements. Also, std (complex number) = sqrt (std(re).^2 + std(im).^2);
	goodstncal = 1;
	stat_std_arr = zeros (1, 6);
	for sta=0:5
		sta_ind = [sta*48+1 : (sta+1)*48];
		% Use only unflagged antennas per station
		stat_std = std(currsol (gainmask(sta_ind) == 0));
		stat_std_re = std(real(currsol (gainmask(sta_ind) == 0)));
		stat_std_im = std(imag(currsol (gainmask(sta_ind) == 0)));

		solwin_std = std (meansol (gainmask(sta_ind) == 0));
		solwin_std_re = std(real (meansol(gainmask(sta_ind) == 0)));
		solwin_std_im = std(imag (meansol(gainmask(sta_ind) == 0)));

		stat_std_arr(sta+1) = stat_std;

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
	if (goodstncal == 1 && err > solparm.errthresh)
		if (mse_ph < 30) % < 30deg. phase err is OK to pass.
			goodstncal = 1;
		else
			goodstncal = 0;
		end;
		% pause;
	end;
	fprintf (1, '{%5.2f, %4.1f, m:%5.2f, r:%5.2f, i:%5.2f, a:%5.2f, p:%5.2f} ', err, pherr, mse, mse_re, mse_im, mse_amp, mse_ph);

	if (debug >= 4) % && goodstncal == 0)
		% Plotting only the bad solutions.
		figure (solparm.gainplt);
		subplot (1, 2, 1);
		% plot (abs(solwindow)');
		plot (real(solwindow)');
		subplot (1, 2, 2);	
		% plot (angle (solwindow)');
		plot (imag (solwindow)');
		
		figure (solparm.currsolplt);
		subplot (1,2,1);
		plot (abs(meansol), '-b');
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
