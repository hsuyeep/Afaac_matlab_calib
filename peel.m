% Script to carry out peeling on a given visibility set. 
% Arguments:
%	acc  : Instantaneous visibility ACM from which to peel
%	tobs : Timestamp of the acc
%  	freq : Frequency stamp of the acc
% srclist: List of sources containing fluxes and positions to peel away.

% Returns:
%	accpeel: The peeled acc.
% pep/23Jan13

function peel (acc, tobs, freq, srclist)
	% For all sources in the source list:

	% Rephase the array towards the source. (To do or not to do?)

	% Generate an Array PSF at the location of the source.

	% Generate model visibilities for the source. 

	% Subtract the model visibilities from original ACM.
