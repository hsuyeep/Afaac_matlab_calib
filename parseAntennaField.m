function [posxpol, posypol, lon, lat, refpos, rotmat] = parseAntennaField(filename, rcumode, EU)

% [posxpol, posypol, lon, lat, refpos, rotmat] =
%     parseAntennaField(filename, rcumode, EU)
%
% Parser for AntennaFields.conf file.
%
% arguments
% filename : filename (including path if necessary) of the AntennaFields.conf
%            file
% rcumode  : RCU mode number (determines arrau configuration)
% EU       : true for European stations
%
% return values
% posxpol : Nelem x 3 matrix with (x, y, z)-positions of x-dipoles (in m)
%           w.r.t. the reference position in ITRF coordinates
% posypol : Nelem x 3 matrix with (x, y, z)-positions of y-dipoles (in m)
%           w.r.t. the reference position in ITRF coordinates
% lon     : geographic longitude of the station (in degrees)
% lat     : geographic latitude of the station (in degrees)
% refpos  : Nelem x 3 vector with (x, y, z)-coordinate of ITRF reference
%           position of the specified antenna array (in m)
% rotmat  : 3 x 3 rotation matrix to convert ITRF to local coordinates
%
% SJW, January 2011
% modified by SJW, January 2012

% parse for antenna locations
if find(rcumode == [1 2 3 4])
    arrayname = 'LBA';
else
    arrayname = 'HBA';
end

fid = fopen(filename, 'r');
token = fscanf(fid, '%s', 1);
prevtoken = [];
while ~((strcmp(token, arrayname) && ~strcmp(prevtoken, 'NORMAL_VECTOR') && ~strcmp(prevtoken, 'ROTATION_MATRIX')) || isempty(token))
    prevtoken = token;
    token = fscanf(fid, '%s', 1);
end
if (strcmp(token, arrayname))
    disp('configuration found, parsing...');
    fscanf(fid, '%s', 2);
    refpos = fscanf(fid, '%f', 3);
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
    return;
end
fclose(fid);

% select the right set of antennas
if EU == 1
    sel = 1:96;
else
    if find(rcumode == [1, 2])
        sel = Nelem/2+1:Nelem;
    else if find(rcumode == [3, 4])
        sel = 1:Nelem/2;
        else
            sel = 1:Nelem;
        end
    end
end
posxpol = posxpol(sel, :);
posypol = posypol(sel, :);

% obtain rotation matrix
fid = fopen(filename, 'r');
token = fscanf(fid, '%s', 1);
prevtoken = [];
while ~((strncmp(token, arrayname, 3) && strcmp(prevtoken, 'ROTATION_MATRIX')) || isempty(token))
    prevtoken = token;
    token = fscanf(fid, '%s', 1);
end
fscanf(fid, '%s', 4);
rotmat = zeros(3, 3);
rotmat(:) = fscanf(fid, '%f', 9);
fclose(fid);

% convert antenna positions to local horizon coordinate system
fid = fopen(filename, 'r');
if strcmp(token, 'HBA0');
    while ~((strcmp(token, 'HBA1') && strcmp(prevtoken, 'ROTATION_MATRIX')) || isempty(token))
        prevtoken = token;
        token = fscanf(fid, '%s', 1);
    end
    fscanf(fid, '%s', 4);
    rotmat2 = zeros(3, 3);
    rotmat2(:) = fscanf(fid, '%f', 9);
    posxpol(1:Nelem/2, :) = posxpol(1:Nelem/2, :) / rotmat;
    posxpol(Nelem/2+1:Nelem, :) = posxpol(Nelem/2+1:Nelem, :) / rotmat2;
    posypol(1:Nelem/2, :) = posypol(1:Nelem/2, :) / rotmat;
    posypol(Nelem/2+1:Nelem, :) = posypol(Nelem/2+1:Nelem, :) / rotmat2;
else
    posxpol = posxpol / rotmat;
    posypol = posypol / rotmat;
end
fclose(fid);

% obtain longitude and latitude of the station
wgs84_f = 1 / 298.257223563;
wgs84_a = 6378137;
wgs84_e2 = wgs84_f * (2 - wgs84_f);
lon = atan2(refpos(2), refpos(1)) * 180 / pi;
r = sqrt(refpos(1)^2 + refpos(2)^2);
prev_lat = 1e5;
lat = atan2(refpos(3), r);
while (abs(lat - prev_lat) >= 1e-12)
    prev_lat = lat;
    normalized_earth_radius = 1 / sqrt((1-wgs84_f)^2 * sin(lat)^2 + cos(lat)^2);
    lat = atan2(wgs84_e2 * wgs84_a * normalized_earth_radius * sin(lat) + refpos(3), r);
end
lat = lat * 180 /pi;