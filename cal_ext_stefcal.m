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
%	sigmas   : Nsrc x 1 vector with initial source power estimates
%	mask     : Nelem x Nelem matrix with ones marking the entries for which
%	           the noise covariance matrix should be estimated and zeros
%	           elsewhere
%   uvflag   : NelemxNelem matrix with ones marking the visibilities to ignore
%			   for other estimators.
%	calim.diffstop: optional argument defining the stop criterion based on the
%	           difference between the solution vectors found in consecutive
%	           iterations, default value is 1e-10
%	calim.maxiter: optional argument defining the maximum number of iterations,
%	           default value is 10
 	
% Return values:
%	g       : Nelem x 1 vector with estimated complex receiver gains
%	sigmas  : Nsrc x 1 vector with estimated apparent source powers
%	Sigma_n : Nelem x Nelem estimated noise covariance matrix
%       --- Returned solution structure; above is obsolete --- 
%	sol.gainsol:Nelem x 1 vector with estimated complex receiver gains from 
%				final iteration.
%	sol.sigmas: Nsrc x 1 vector with estimated apparent source powers from 
%			    final iteration.
%	sol.sigman: Nelem X Nelem estimated noise covariance matrix.
%	sol.calext_iters: Number of cal_ext iterations for convergence.
%	sol.pinv_sol: Convergence condition on final iteration.
%   sol_gainsolv: Solution structure from the gainsolv call in the *final* 
%				  cal_ext. iteration.
 
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
%
%  Modified to return a structured gain solution, including stefcal gains.

% function [g, sigmas, Sigma_n] = cal_ext_stefcal(Rhat, A, sigmas, mask, calim)
function [sol, sol_gainsolv] = cal_ext_stefcal(Rhat, A, sigmas, mask, uvflag, calim, ptsrc, up)
	% parameters
	Nelem = size(A,1);
    Nsrc  = length (sigmas);
    Ndir  = size (calim.ncomp2src);
	
	% initialization
	ghat = zeros(Nelem, calim.maxiter+1);    % Scalar gain
	ghat(:, 1) = 1;                          % Set first iteration to unit gain.
	sigmahat = zeros(Nsrc, calim.maxiter+1); % Model src. flux estimates.

    % Incoming model src. flux estimates, set invalid entries to 1.
	sigmas(~isfinite(sigmas)) = 1;           

	sigmahat(:, 1) = sigmas;                 % Initialize model src. estimate.
	Sigma_n = zeros(Nelem);					 % Initialize noise matrix.
	
	total_stefcal_iters = 0;
	vismask = 1-(mask | uvflag);             % 1's in both mask & uvflag mean 
											 % ignore the baseline.
	% implementation using WALS
	for iter = 2:calim.maxiter+1
	    
		%------------ Per sensor gain estimation -------------
	    % estimate g using baseline restriction (eq. 36 of [1])
		%--------- UNComment this part to get to the full WALS solver --------%
%	    alpha = computeAlphaW(Rhat .* (1 - mask), ... 
%					(A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask));
%	    if (sum(~isfinite(alpha(:))) == 0)
%	       [v, d] = eig(alpha);
%	       [~, idx] = max(diag(d));
%	       ghat(:, iter) = conj(v(:, idx));
%	    else
%	       ghat(:, iter) = 1;
%	    end
%		sol_gainsolv.ghat = ghat(:, iter);	
%		sol_gainsolv.gresult = 1;
%		sol_gainsolv.iter = 1;
%		sol_gainsolv.tol = 1;
		%--------- UNComment this part to get to the full WALS solver --------%

	    
		%------------ Per sensor gain estimation -------------
	    % gain estimation using StefCal - experimental
		%--------- Comment this part to get to the full WALS solver --------%

        % Assuming same number of components per model source.
        % srcflux = reshape (repmat (sigmahat(:, iter-1), [1,Ndir])', [1,Ndir*Nsrc]);

        srcflux = [];
        if (ptsrc == 0)
            for src = 1:length (up)
                if (up(src) == 1)
                    valid_components = calim.val_comp{src};
                    apparentmodelflux = sigmahat(src, iter-1) .* calim.modnormimg{src}(valid_components);
                    srcflux = [srcflux; apparentmodelflux];
                end;
            end;
        else
            srcflux = sigmahat (:, iter-1);
        end;

        Rmodel = (A * diag(srcflux) * A'); 
 	    sol_gainsolv = gainsolv (1e-6, Rmodel ...
 			.* vismask, Rhat .* vismask, ghat(:, iter-1),...  
			calim.maxiter_gainsolv, calim.debug);
		total_stefcal_iters = total_stefcal_iters + sol_gainsolv.iter;
 	    ghat(:, iter) = sol_gainsolv.ghat;
		%--------- Comment this part to get to the full WALS solver --------%

%	    GA = diag(ghat(:, iter)) * A;
%	    Rest = GA * diag(sigmahat(:, iter-1)) * GA';
%	    compidx = (1 - mask) ~= 0;
%
%	    normg = abs(sqrt(pinv(Rest(compidx)) * Rhat(compidx)));
%	    ghat(:,iter) = normg * ghat(:, iter) / (ghat(1, iter)  ... 
%					   / abs(ghat(1, iter)));
%	    ghat(~isfinite(ghat(:, iter)), iter) = 1;

		% NEW CONSTRAINT: Avg. of all gain solutions is unity.
		avg_gain = mean (abs(ghat(:, iter)));
		ghat (:, iter) = ghat (:,iter)/avg_gain;
		% Phase constraint: Phase of first antenna is 0.
		ghat (:, iter) = ghat (:,iter) / (ghat (1, iter)/abs(ghat(1, iter)));
	
		%------------ Model source flux estimation -------------
	    % estimate sigmahat using sigmanhat and ghat (eq. 42 of [1])
        % fprintf (1, 'cal_ext: iter=%d, Rank/rcond(Rhat)=%f/%f\n', ...
		% 			iter, rank(Rhat), rcond(Rhat));
	    invR = inv(Rhat);
        
        diravgA = [];
        if (ptsrc == 0)
            cnt = 1;
            for src = 1:Nsrc
                diravgA = [diravgA; mean (A(:,cnt:calim.ncomp2src(src), 2))];
                cnt = cnt + calim.ncomp2src(src);
            end;
        else
            diravgA = A;
        end;
        
        % Used if the number of directions is exactly the same for all sources.
        % diravgA = squeeze(mean (reshape (A, [Nelem, Ndir, Nsrc]), 2));
	    GA = diag(ghat(:, iter)) * diravgA; % new line, use normalized G
%	    sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) ... 
%  						* diag(GA' * invR * (Rhat - Sigma_n) * invR * GA));
%%  						* diag(GA' * invR * (Rhat .* (1-mask)) * invR * GA));
		r = Rhat (mask(:) == 0);   % Ignore the shortest baselines, take the others.
		M = khatrirao(conj(GA), GA);
		M = M(mask(:)==0, :);
		sigmahat(:,iter) = real ((M'*M)\M' * r);

	    if sum(~isfinite(sigmahat(:, iter))) ~= 0
	        sigmahat(:, iter) = sigmahat(:, iter-1);
	    end
	
	    % use first source as amplitude reference (normalize by CasA flux).
%	    if (sigmahat(1, iter) ~= 0)
%	        sigmahat(:, iter) = sigmahat(:, iter) / sigmahat(1, iter);
%	    end
	
	    % remove negative values
%	    sigmahat(:, iter) = max(sigmahat(:, iter), 0);
	    
		%------------ Noise covariance estimation -------------
		%i----------- assuming Sigma_n = diag(sig_n)-------------
	    % estimate sigmanhat using sigmahat and ghat (eq. 44 of [1])
	    Sigma_n = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* mask;
	    % Sigma_n = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* (1-mask);
	
	    % test for convergence
	    theta_prev = [ghat(:, iter-1); sigmahat(:, iter-1)];
	    theta = [ghat(:, iter); sigmahat(:, iter)];
		pinv_sol = abs(pinv(theta_prev) * theta - 1);
	    % if (abs(pinv(theta_prev) * theta - 1) < calim.diffstop)

		if (calim.debug > 1)
			fprintf (1, 'cal_ext: Iteration %d\n', iter);
		end;
	    if (pinv_sol < calim.diffstop)
	        break;
	    end
	end
	
	if (calim.debug > 0)
	  fprintf(1,'Cal_ext_stefcal:stop after %d iter, pinv: %.4f.\n', ... 
				iter, pinv_sol);
	end
	sol.gainsol = (1 ./ ghat(:, iter))';   % NOTE: inverse of estimated gain!
	sol.sigmas = sigmahat(:, iter);
	sol.sigman = Sigma_n;
	sol.calext_iters = iter;
	sol.pinv_sol = pinv_sol;
	sol.stefcal_iters = total_stefcal_iters;
	sol.stefcal_tol = sol_gainsolv.tol;
