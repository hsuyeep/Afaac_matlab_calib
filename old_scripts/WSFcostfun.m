function cost = WSFcostfun(theta, EsWEs, G, freq, xpos, ypos)

% cost = WSFcostfun(theta, EsWs, G, freq, xpos, ypos)
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
% xpos/ypos: Nelem x 1 matrix with ITRF positions of array elements, rotated
%           to make CS002 as the reference. zpos is unused due to a coplanar
%           array.
%
% Return value
% cost : value of the cost function at theta
%
% SJW, Unknown date, comment added by Peeyush 24Jul12
nsrc = length(theta)/2;
c = 2.99792e8;
nelem = length(xpos);
lsrc = theta(1:nsrc);
msrc = theta(nsrc+1:end);

A = G * exp(-2 * pi * 1i * freq/c * (xpos.' * lsrc.' + ypos.' * msrc.'));
PAperp = eye(nelem) - A * inv(A' * A) * A';
cost = trace(PAperp * EsWEs);
