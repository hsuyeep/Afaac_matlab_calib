% Script to carry out tracking  calibration. 
% Modification of cal_ext_stefcal.m. Main difference is
% carrying out a single iteration of calibration, using previous solutions as
% a starting point.
% pep/25Jan13.

% Array calibration assuming the model
% R = G * A * Sigma_s * A' * G' + Sigma_n
% where
% G       : diagonal Nelem x Nelem gain matrix with complex receiver gains
%           on the main diagonal, assumed unknown.
% A       : Nelem x Nsrc matrix describing the geometrical delays, assumed
%           known.
% Sigma_s : diagonal Nsrc x Nsrc matrix with apparent source powers on the
%           main diagonal, assumed unknown.
% Sigma_n : Nelem x Nelem noise covariance matrix, which may have a number
%           of non-zero (off-diagonal) elements, which are assumed unknown.
%
% The calibration is done using the weighted alternating least squares
% method described in [1] with extension to non-zero noise covariance
% matrix as described in [2].
%
% Arguments:
%	Rhat     : Nelem x Nelem measured array covariance matrix
%	A        : Nelem x Nsrc known geometrical delays
%   sel		 : Vector containing mask for sky-model sources actually used.
%	mask     : Nelem x Nelem matrix with ones marking the entries for which
%	           the noise covariance matrix should be estimated and zeros
%	           elsewhere
%	calim.diffstop: optional argument defining the stop criterion based on the
%	           difference between the solution vectors found in consecutive
%	           iterations, default value is 1e-10
%	calim.maxiter: optional argument defining the maximum number of iterations,
%	           default value is 10
%	prevsol : Structure containing best estimates from a previous timeslice.
 	
% Return values:
%	sol.g       : Nelem x 1 vector with estimated complex receiver gains
%	sol.sigmas  : size (rodata.srclist) x 1 vector with estimated apparent 
%				  source powers, at the location of selected sources.
%	sol.Sigma_n : Nelem x Nelem estimated noise covariance matrix
 
% References:
% 	[1] Stefan J. Wijnholds and Alle-Jan van der Veen, "Multisource
% 		Self-calibration for Sensor Arrays", IEEE Transactions on Signal
% 		Processing, V57, no. 9, pp3512-3522, September 2009
% 	[2] Stefan J. wijnholds and Alle-Jan van der Veen, "Self-Calibration of
% 		Radio Astronomical Arrays With Non-Diagonal Noise Covariance Matrix",
% 		17th European Signal Processing Conference, August 24-28, 2009, Glasgow
% 		(UK)
%
% SJW, 16 June 2010
% Routine modified to utilize the StefCal implementation provided by SJW.
% Pep, sometime Mar, 2012

function [sol, sol_gainsolv] = track_cal_ext_stefcal(Rhat, A, sel, mask, ... 
								calim, prevsol)
	% parameters
	[Nelem, Nsrc] = size(A);
	
	% initialization now from prevsol.
	ghat = zeros(Nelem, calim.maxiter+1);    % Scalar gain

	% Set first iteration to prev. estimated gain.
 	ghat(:, 1) = 1 ./ prevsol.gainsol'; 
	% ghat(:, 1) = 1;                        % Set first iteration to unit gain.
	sigmahat = zeros(Nsrc, calim.maxiter+1); % Model src. flux estimates.
	sigmahat(:, 1) = prevsol.sigmas(sel);    % Initialize model src. estimate.

	sol.sigman = prevsol.sigman; %zeros(Nelem);	   % Initialize noise matrix.
	
	if (calim.debug > 0)
		fprintf (1, 'track_cal_ext: Sigmas from prev. solution: \n');
		disp (sigmahat(:,1)');
	end;

    % Incoming model src. flux estimates, set invalid entries to 1.
	% sigmas(~isfinite(sigmas)) = 1;           

% % TO EXPT. WITH CONVERGENCE, BUT WITH INITIAL GUESSES FROM PREV. SOLU.
%	ghat = zeros(Nelem, 2);    % Scalar gain, current and previous values.
% 	ghat(:, 1) = 1 ./ prevsol.gainsol; % Set first iteration to prev. estimated gain.
% 	% ghat(:, 1) = 1; % Set first iteration gains to 1.
%
%	sigmahat = zeros(Nsrc, 2); % Model src. flux estimates.
%	sigmahat(:, 1) = prevsol.sigmas(sel);   % Initialize model src. estimate.
%	% Sigma_n = zeros(Nelem);	   % Initialize noise matrix.
%	Sigma_n = prevsol.sigman;	   % Initialize noise matrix with prev. sol.
%
	
	% implementation of a single iteration of calibration using WALS
	total_stefcal_iters = 0;
	iter = 2;
	% for iter = 2:calim.maxiter+1
	    
		%------------ Per sensor gain estimation -------------
	    % estimate g using baseline restriction (eq. 36 of [1])
	    % alpha = computeAlphaW(Rhat .* (1 - mask), ... 
		%			(A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask));
	    % if (sum(~isfinite(alpha(:))) == 0)
	    %    [v, d] = eig(alpha);
	    %    [~, idx] = max(diag(d));
	    %    ghat(:, iter) = conj(v(:, idx));
	    % else
	    %    ghat(:, iter) = 1;
	    % end
	    
		%------------ Per sensor gain estimation -------------
	    % gain estimation using StefCal - experimental
 	    sol_gainsolv = gainsolv (1e-6, (A * diag(sigmahat(:, iter-1)) * A')...
 						.* (1 - mask), Rhat .* (1 - mask), ghat(:, iter-1), ...
						calim.maxiter_gainsolv, calim.debug);
		total_stefcal_iters = total_stefcal_iters + sol_gainsolv.iter;
		ghat(:, iter) = sol_gainsolv.ghat;

	    GA = diag(ghat(:, iter)) * A;
	    Rest = GA * diag(sigmahat(:, iter-1)) * GA';
	    compidx = (1 - mask) ~= 0;

	    normg = abs(sqrt(pinv(Rest(compidx)) * Rhat(compidx)));
	    ghat(:,iter) = normg * ghat(:, iter) / (ghat(1, iter)  ... 
					   / abs(ghat(1, iter)));
	    ghat(~isfinite(ghat(:, iter)), iter) = 1;
	
		%------------ Model source flux estimation -------------
	    % estimate sigmahat using sigmanhat and ghat (eq. 42 of [1])
        % disp (['cal_ext: iter= ' num2str(iter) ', Rank/rcond(Rhat) = ' num2str(rank(Rhat)) ' ' num2str(rcond(Rhat))]);
	    invR = inv(Rhat);
        
	    GA = diag(ghat(:, iter)) * A; % new line, use normalized G
	    sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) ... 
							* diag(GA' * invR * (Rhat - sol.sigman) * invR * GA));
	    if sum(~isfinite(sigmahat(:, iter))) ~= 0
	        sigmahat(:, iter) = sigmahat(:, iter-1);
	    end
	
	    % use first source as amplitude reference (normalize by CasA flux).
	    if (sigmahat(1, iter) ~= 0)
	        sigmahat(:, iter) = sigmahat(:, iter) / sigmahat(1, iter);
	    end
	
	    % remove negative values
	    sigmahat(:, iter) = max(sigmahat(:, iter), 0);
	    
		%------------ Noise covariance estimation -------------
		%i----------- assuming Sigma_n = diag(sig_n)-------------
	    % estimate sigmanhat using sigmahat and ghat (eq. 44 of [1])
	    sol.sigman = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* mask;
	
	    % test for convergence
	    theta_prev = [ghat(:, iter-1); sigmahat(:, iter-1)];
	    theta = [ghat(:, iter); sigmahat(:, iter)];
		pinv_sol = abs(pinv(theta_prev) * theta - 1);
%	    if (pinv_sol < calim.diffstop)
%	        break;
%	    end
	    fprintf (1, 'track_cal_ext: Iteration difference: %g.\n', pinv_sol);
 	% end
	
	sol.pinv_sol = pinv_sol;
	sol.calext_iters = calim.maxiter;
	sol.gainsol = (1 ./ ghat(:, iter))';
	sol.sigmas = zeros (size (prevsol.sigmas));
	sol.sigmas(sel) = sigmahat(:, iter);
	sol.stefcal_iters = total_stefcal_iters;
	sol.stefcal_tol = sol_gainsolv.tol;
