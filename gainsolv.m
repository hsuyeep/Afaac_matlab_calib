%
% Gain estimation using algorithm developed by Stef Salvini, edited by SJW
% for readability.
%
% Arguments:
%	tolopt : stop criterion on the change in the norm of the estimated gain
%	         vectors in consecutive iterations
%	R0     : Nelem x Nelem array covariance matrix mode, zero entries are
%	         flagged
%	Rhat   : Nelem x Nelem measured array covariance matrix, zero entries are
%	         flagged
%	g0     : Nelem x 1 complex vector with initial estimate for gains
% maxiter  : Maximum number of iterations to carry out. Loop breaks on either 
%			 the maxiter or tolopt. option being satisfied, whichever occurs 
%			 first.
%   debug  : Debug level for gainsolv.
%
% Return values:
%	sol.ghat   : Final direction independent gain solutions, NelemX1.
%	sol.gresult: History of gain solutions over iterations.
%	sol.iter   : Number of iterations
%	sol.tol    : Final achieved tolerance.
%   sol.dgnorm : gain norm derivative as a function of iteration number.
%   sol.gnorm  : gain norm as a function of iteration number.
%
% SS & SJW, March 7, 2012
% Per iteration solution extraction by pep, 05Feb13 

% function g = gainsolv(tolopt, R0, Rhat, g0)
function [sol] = gainsolv(tolopt, R0, Rhat, g0, maxiter, debug)

	% initialization
	Nelem = length(g0);
	g = g0;
	colnorm = zeros(Nelem, 1);
	Rhatnorm = Rhat;             % Initialize with measured ACM.
	dgnorm = 0;
	gnorm = 0;

	% normalize the measured array covariance matrix (weighting)
	for idx = 1:Nelem
	    colnorm(idx) = Rhat(:, idx)' * Rhat(:, idx);
	    Rhatnorm(:, idx) = Rhat(:, idx) / colnorm(idx);
	end

	% initial calibration
	Rcolcal = complex(zeros(Nelem), zeros(Nelem));
	for idx = 1:Nelem
	    Rcolcal(:, idx) = g .* R0(:, idx);
	end

%	fprintf (2, 'gainsolv: gnorm: %f\n', norm (g));
%	figure;
%	subplot (1,2,1);
%	plot (abs(g));
%	subplot (1,2,2);
%	plot (angle(g));
%	drawnow();
%	pause;

	% start iterating
	calest = complex(zeros(Nelem, 1), zeros(Nelem, 1));
	gresult = zeros(Nelem, maxiter);
	sol.gnorm = zeros (1, maxiter);
	sol.dgnorm = zeros (1, maxiter);

	for iter = 1:maxiter
	    % update gain estimate
	    for idx = 1:Nelem
	        calest(idx) = Rhatnorm(:, idx)' * Rcolcal(:, idx);
	    end
	    g = 1 ./ conj(calest);
        %disp(['iter: ' num2str(iter) ' Rank of calest:' num2str(rank(calest))]);
	    if (mod(iter, 2) > 0)
	        % find new convergence curve
	        gest1 = g;
	        % dgnorm = norm(g - gold);
	        % gnorm = norm(g);
	    else
	        % update calibration results so far and check for convergence
	        gold = g;
	        g = (g + gest1) / 2;
	        if (iter >= 2)
	            dgnorm = norm(g - gold);
	            gnorm = norm(g);
				if (debug > 1)
					fprintf (1,'gainsolv: 	iter %2d: dgnorm: %.6f, gnorm: %.4f.\n',...
							 iter, dgnorm, gnorm);
				end;
	            if (dgnorm / gnorm <= tolopt)
	                break;
	            end
	        end
	    end
		
		% Update statistics on every iteration
		sol.gnorm (iter)  = gnorm;
		sol.dgnorm (iter) = dgnorm;

	    % update calibration of array covariance matrix
	    for idx = 1:Nelem
	        Rcolcal(:, idx) = g .* R0(:, idx);
	    end
	    gresult(:, iter) = g;
	end
	if (debug > 0)
		fprintf(1,'gainsolv: iters: %2d, dgnorm: %.4f, gnorm: %.4f.\n\n',...
					iter, dgnorm, gnorm);
	end;
	sol.ghat = g;           % Final direction independent gain solution.
	sol.gresult = gresult;  % History of gain solutions over iterations.
	sol.iter = iter;        % Number of iterations
	sol.tol = dgnorm/gnorm; % Final achieved tolerance.
