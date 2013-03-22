function XYPhaseDelta(calinfile, XYPhase, caloutfile)

% XYPhaseDelta(calinfile, XYPhase, caloutfile)
%
% Change the phase offset between the array of x-dipoles and the array of
% y-dipoles in the calibration table. It is assumed by convention that the
% x-dipoles are connected to even numbered RCUs (when starting from 0) and
% the y-dipoles are connected to add numbered RCUs.
%
% This reads the calibration table, modifies the phases and writes the
% calibration table back to disk. If calinfile and caloutfile are the same,
% the calibration table file will be overwritten!
%
% Arguments
% calinfile  : name of the file with the calibration table to be corrected
% XYPhase    : phase change to be applied (in radians)
% caloutfile : name of the file in which to write the new calibration table
%
% SJW, 18 October 2011

% read calibration table
[calx, caly] = readCalTable(calinfile);

% apply phase correction
caly = caly * exp(1i * XYPhase);

% write new calibration table
writeCalTable(calx, caly, caloutfile);
