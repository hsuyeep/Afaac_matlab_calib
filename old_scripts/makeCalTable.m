% makeCaltable
%
% This script faciltiates the construction of a calibration table based on
% the results produced by the calibrateData routine. The latter produces a
% struct called "calres" and it is assumed that this struct is loaded to
% the Matlab workspace.
%
% SJW, December 2010
% modified by SJW, October 2011

%% start reduction of x-polarization

% determine the number of scans in the calibration observation (Nobs), the
% number of subbands (Nfreq) and the number of antennas (Nant)
[Nobs, Nfreq, Nant] = size(calres.calx);

% specify polarization fsor which to produce the calibration table
cal = calres.calx;
 
% collect center frequency of each subband
% the 10th scan is used to avoid problems with slow starting stations
freq = calres.freq(2, :);

% initially, we assume that no subbands need to be flagged manually
sbflag{Nant} = [];

% fit amplitude model
fitAmplModel;

% fit phase model
fitPhaseModel;

%% collect result for x-polarization in workspace variable
% combine amplitude and phase results and fill in uncalibrated entries with
% a reasonable value
calx = diag(abscal) * exp(1i * argcal);
calx(isnan(calx)) = 1;

%% repeat the exercise for y-polarization
cal = calres.caly;
sbflag{Nant} = [];
fitAmplModel;
fitPhaseModel;

%% collect result for y-polarization in workspace variable
caly = diag(abscal) * exp(1i * argcal);
caly(isnan(caly)) = 1;
