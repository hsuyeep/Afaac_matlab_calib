function [posxpol, posypol, lon, lat, height] = parseAntennaArrays(filename, arrayname)

% [posxpol, posypol, lon, lat, height] = 
%     parseAntennaArrays(filename, arrayname)
%
% Parser for AntennaArrays.conf file
%
% arguments
% filename  : filename (including path if necessary) of the
%             AntennaArrays.conf file
% arrayname : name of the array for which information has to be extracted
%             from the configuration file
%
% return values
% posxpol : Nelem x 3 matrix with the (x, y, z)-positions of x-dipoles
% posypol : Nelem x 3 matrix with the (x, y, z)-positions of y-dipoles
% lon     : geographic longitude of the station
% lat     : geographic latitude of the station
% height  : geographic altitude of the station
%
% SJW, 2009

fid = fopen(filename, 'r');
token = fscanf(fid, '%s', 1);
while ~(strcmp(token, arrayname) || isempty(token))
    token = fscanf(fid, '%s', 1);
end
if (strcmp(token, arrayname))
    disp('configuration found, parsing...');
    % obtain station location
    fscanf(fid, '%s', 2);
    lon = fscanf(fid, '%f', 1);
    lat = fscanf(fid, '%f', 1);
    height = fscanf(fid, '%f', 1);
    fscanf(fid, '%s', 1);
 
    % obtain antenna locations
    Nelem = str2double(fscanf(fid, '%s', 1));
    fscanf(fid, '%s', 5);
    posxpol = zeros(Nelem, 3);
    posypol = zeros(Nelem, 3);
    for idx = 1:Nelem
        posxpol(idx, :) = fscanf(fid, '%f', 3);
        posypol(idx, :) = fscanf(fid, '%f', 3);
    end
else
    disp('configuration not found');
    posxpol = [];
    posypol = [];
    lon = [];
    lat = [];
    height = [];
end
