function [posx, posy, normal] = parseAntennaFieldITRF(filename, rcumode, EU)

% [posx, posy, normal] = parseAntennaFieldITRF(filename, rcumode, EU)
%
% Parse antenna positions in ITRF and the normal vector to the station
% field from the AntennaFields.conf file.
%
% arguments
% filename : filename (including path if necessary) of the
%            AntennaFields.conf file
% rcumode  : RCU mode number (determines array configuration)
% EU       : true for European stations
%
% return values
% posx   : Nelem x 3 matrix with (x, y, z)-positions of x-dipoles (in m)
% posy   : Nelem x 3 matrix with (x, y, z)-positions of y-dipoles (in m)
% normal : 3 x 1 normal vector to the station field
%
% SJW, 18 May 2011

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
    disp('configuration found, parsing ...');
    fscanf(fid, '%s', 6);
    % obtain antenna locations
    Nelem = str2double(fscanf(fid, '%s', 1));
    fscanf(fid, '%s', 5);
    posx = zeros(Nelem, 3);
    posy = zeros(Nelem, 3);
    for idx = 1:Nelem
        posx(idx, :) = fscanf(fid, '%f', 3);
        posy(idx, :) = fscanf(fid, '%f', 3);
    end
else
    disp('configuration not found');
    posx = [];
    posy = [];
    normal = [];
    return;
end
fclose(fid);

% select the right set of antennas
if EU == 1
    sel = 1:96;
else
    if find(rcumode == [1 2])
        sel = Nelem/2+1:Nelem;
    else if find(rcumode == [3, 4])
            sel = 1:Nelem/2;
        else
            sel = 1:Nelem;
        end
    end
end
posx = posx(sel, :);
posy = posy(sel, :);

% obtain normal vector
fid = fopen(filename, 'r');
token = fscanf(fid, '%s', 1);
prevtoken = [];
while ~((strncmp(token, arrayname, 3) && strcmp(prevtoken, 'NORMAL_VECTOR')) || isempty(token))
    prevtoken = token;
    token = fscanf(fid, '%s', 1);
end
fscanf(fid, '%s', 2);
normal = fscanf(fid, '%f', 3);
