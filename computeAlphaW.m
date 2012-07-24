function alphamat = computeAlphaW(Rhat, R0)

% alphamat = computeAlpha(Rhat, R0)
%
% Support function for cal_msrc which computes a matrix
% alphamat = g * (1./g).' which corresponds to the matrix M in [1] and [2], 
% where g are the complex gains of the elements producing the input signals for
% which Rhat is the measured array covariance matrix (ACM) and R0 is the
% model ACM. If specific entries in these matrices are set to zero, these
% entries are not used for estimating alpha. This feature can be exploited
% when a baseline restriction is applied.
%
% Parameters 
% Rhat     : P x P measured ACM
% R0       : P x P expected (model) ACM
%
% Return value
% alphamat : P x P matrix representing the best estimate for g * (1./g).'
%
% References
% [1] Stefan J. Wijnholds and Albert-Jan Boonstra, "A Multisource
% Calibration Method for Phased Array Radio Telescopes", IEEE Workshop on
% Sensor Array and Multichannel Processing (SAM), Waltham (MA), USA,
% July 2006
% [2] Stefan J. Wijnholds and Alle-Jan van der Veen, "Multisource Self-
% Calibration for Sensor Arrays", IEEE Transactions on Signal Processing, V57,
% no. 9, pp3512-3522, September 2009
%
% SJW, 2006

% determine the number of elements
[~, nelem] = size(Rhat);
% ignore the autocorrelations since these are the sum of source and
% receiver powers
R0_red = R0 - diag(diag(R0));

alphamat = zeros(nelem, nelem);
for i = 1:nelem
    for j = 1:nelem
        nonzero = (R0_red(i, :) ~= 0) & (R0_red(j, :) ~= 0) & (Rhat(i, :) ~= 0) & (Rhat(j, :) ~= 0);
        rhati = Rhat(i, nonzero)';
        rhatj = Rhat(j, nonzero)';
        r0i = R0(i, nonzero)';
        r0j = R0(j, nonzero)';
        w = abs(r0i .* rhatj).^2;
        alphamat(i, j) =  sum(w .* (rhati .* r0j) ./ (rhatj .* r0i)) / sum(w);
    end
end
