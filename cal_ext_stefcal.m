function [g, sigmas, Sigma_n] = cal_ext_stefcal(Rhat, A, sigmas, mask, debug,  varargin)

% [g, sigmas, Sigma_n] = cal_ext(Rhat, A, sigmas, mask, diffstop, maxiter)
%
% Array calibration assuming the model
%
% R = G * A * Sigma_s * A' * G' + Sigma_n
%
% where
% G       : diagonal Nelem x Nelem gain matrix with complex receiver gains
%           on the main diagonal, assumed unknwon
% A       : Nelem x Nsrc matrix describing the geometrical delays, assumed
%           known
% Sigma_s : diagonal Nsrc x Nsrc matrix with apparent source powers on the
%           main diagonal, assumed unknown
% Sigma_n : Nelem x Nelem noise covariance matrix, which may have a number
%           of non-zero (off-diagonal) elements, which are assumed unknown
%
% The calibration is done using the weighted alternating least squares
% method described in [1] with extension to non-zero noise covariance
% matrix as described in [2].
%
% arguments
% Rhat     : Nelem x Nelem measured array covariance matrix
% A        : Nelem x Nsrc known geometrical delays
% sigmas   : Nsrc x 1 vector with initial source power estimates
% mask     : Nelem x Nelem matrix with ones marking the entries for which
%            the noise covariance matrix should be estimated and zeros
%            elsewhere
% diffstop : optional argument defining the stop criterion based on the
%            difference between the solution vectors found in consecutive
%            iterations, default value is 1e-10
% maxiter  : optional argument defining the maximum number of iterations,
%            default value is 10
%
% return values
% g       : Nelem x 1 vector with estimated complex receiver gains
% sigmas  : Nsrc x 1 vector with estimated apparent source powers
% Sigma_n : Nelem x Nelem estimated noise covariance matrix
%
% References
% [1] Stefan J. Wijnholds and Alle-Jan van der Veen, "Multisource
% Self-calibration for Sensor Arrays", IEEE Transactions on Signal
% Processing, V57, no. 9, pp3512-3522, September 2009
% [2] Stefan J. wijnholds and Alle-Jan van der Veen, "Self-Calibration of
% Radio Astronomical Arrays With Non-Diagonal Noise Covariance Matrix",
% 17th European Signal Processing Conference, August 24-28, 2009, Glasgow
% (UK)
%
% SJW, 16 June 2010
% Routine modified to utilize the StefCal implementation provided by SJW.
% Pep, sometime Mar, 2012

% parameters
[Nelem, Nsrc] = size(A);
if (nargin > 5)
    diffstop = varargin{1};
else
    diffstop = 1e-10;
end
if (nargin > 6)
    maxiter = varargin{2};
else
    maxiter = 10;
end

% initialization
ghat = zeros(Nelem, maxiter+1);
ghat(:, 1) = 1;
sigmahat = zeros(Nsrc, maxiter+1);
sigmas(~isfinite(sigmas)) = 1;
sigmahat(:, 1) = sigmas;
Sigma_n = zeros(Nelem);

% implementation using WALS
for iter = 2:maxiter+1
    
    % original code
    % estimate g using baseline restriction
    %alpha = computeAlphaW(Rhat .* (1 - mask), (A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask));
    %if (sum(~isfinite(alpha(:))) == 0)
    %    [v, d] = eig(alpha);
    %    [~, idx] = max(diag(d));
    %    ghat(:, iter) = conj(v(:, idx));
    %else
    %    ghat(:, iter) = 1;
    %end
    
    % gain estimation using StefCal - experimental
    ghat(:, iter) = gainsolv(1e-6, (A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask), Rhat .* (1 - mask), ghat(:, iter-1));
    GA = diag(ghat(:, iter)) * A;
    Rest = GA * diag(sigmahat(:, iter-1)) * GA';
    compidx = (1 - mask) ~= 0;
    normg = abs(sqrt(pinv(Rest(compidx)) * Rhat(compidx)));
    ghat(:, iter) = normg * ghat(:, iter) / (ghat(1, iter) / abs(ghat(1, iter)));
    ghat(~isfinite(ghat(:, iter)), iter) = 1;

    % estimate sigmahat using sigmanhat and ghat
    invR = inv(Rhat);
    GA = diag(ghat(:, iter)) * A; % new line, use normalized G
    sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) * diag(GA' * invR * (Rhat - Sigma_n) * invR * GA));
    if sum(~isfinite(sigmahat(:, iter))) ~= 0
        sigmahat(:, iter) = sigmahat(:, iter-1);
    end
    % use first source as amplitude reference
    if (sigmahat(1, iter) ~= 0)
        sigmahat(:, iter) = sigmahat(:, iter) / sigmahat(1, iter);
    end
    % remove negative values
    sigmahat(:, iter) = max(sigmahat(:, iter), 0);
    
    % estimate sigmanhat using sigmahat and ghat
    Sigma_n = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* mask;

    % test for convergence
    theta_prev = [ghat(:, iter-1); sigmahat(:, iter-1)];
    theta = [ghat(:, iter); sigmahat(:, iter)];
    if (abs(pinv(theta_prev) * theta - 1) < diffstop)
        break
    end
end
if (debug > 0)
  disp(['Cal_ext_stefcal:optimization stopped after ' num2str(iter) ' iteration(s).']);
end
g = ghat(:, iter);
sigmas = sigmahat(:, iter);
