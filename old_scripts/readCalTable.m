function [calx, caly] = readCalTable(calfile)

% [calx, caly] = readCalTable(calfile)
%
% Read calibration table from a properly formatted calibration table file.
%
% Argument
% calfile : file name of calibration table file
%
% Return values
% calx : Nelem x Nsb (48 x 512 for Dutch stations, 96 x 512 for European
%        stations) matrix containing a complex valued correction factor for
%        each element and each subband for the array of x-elements.
% caly : corresponding result for the array of y-elements.
%
% SJW, 18 October 2011

% read calibration table
fid = fopen(calfile, 'r');
data = fread(fid, 'double');
fclose(fid);
Nrcu = length(data) / (2 * 512);
cal = reshape(data(1:2:end) + 1i * data(2:2:end), [Nrcu, 512]);
calx = cal(1:2:end, :);
caly = cal(2:2:end, :);
