function cost = WSFcostfunITRF(theta, EsWEs, G, freq, posITRF)

% cost = WSFcostfunITRF(theta, EsWs, G, freq, posITRF)
%
% Compute the value of the cost function for Weighted Subspace Fitting.
%
% Arguments
% theta   : 2 * Nsrc-element vector with parameter values for which to
%           evaluate the cost function. The first Nsrc elements represent
%           azimuthal angles in radian, the second Nsrc elements are
%           interpreted as elevation angles in radian.
% EsWEs   : Nelem x Nelem matrix formed by the weighted signal subspace
% G       : Nelem x Nelem diagonal gain matrix with element gains
% freq    : frequency in Hz
% posITRF : Nelem x 3 matrix with ITRF positions of array elements
%
% Return value
% cost : value of the cost function at theta
%
% SJW, 27 March 2012

nsrc = length(theta)/2;
c = 2.99792e8;
nelem = size(posITRF, 1);
phisrc = theta(1:nsrc);
thsrc = theta(nsrc+1:end);
srcpos = [cos(phisrc) .* cos(thsrc), sin(phisrc) .* cos(thsrc), sin(thsrc)];

A = G * exp(-2 * pi * 1i * freq/c * (posITRF * srcpos.'));
PAperp = eye(nelem) - A * inv(A' * A) * A';
cost = trace(PAperp * EsWEs);
