function [calx, caly, header] = readCalTable(calfile)

% [calx, caly, header] = readCalTable(calfile)
%
% Read calibration table from a properly formatted calibration table file.
%
% Argument
% calfile : file name of calibration table file
%
% Return values
% calx   : Nelem x Nsb (48 x 512 for Dutch stations, 96 x 512 for European
%          stations) matrix containing a complex valued correction factor
%          for each element and each subband for the array of x-elements.
% caly   : corresponding result for the array of y-elements.
% header : header of calibration table
%
% SJW, 18 October 2011
% modified by SJW on 3 February 2012 (header in calibration tables)
% modified by ME on 27 February 2012 (according to latest header
% definitions.
% N.B. this version is able to read CalTables with or with no headers.

% get file descriptor
fid          = fopen(calfile, 'r');

% distinguish between CalTable files with or without headers; newer
% CalTable files should contain headers
header       = [];
header_start = 'HeaderStart';
header_stop  = 'HeaderStop';

newline      = fgetl(fid);

if strcmp(newline,header_start)
    position_indicator     = ftell(fid);
    newline                = fgetl(fid);
    fseek(fid,position_indicator,-1);
    
    while (~strcmp(newline, header_stop))    
        curchar = fscanf(fid, '%c', 1);
        header = [header, curchar];
        
        if int8(curchar) == 10                                             % 10 is integer value of character '\n' (new line)
            position_indicator = ftell(fid);
            newline            = fgetl(fid);
            fseek(fid,position_indicator,-1);
        end
    end
    newline                = fgetl(fid);                                   % dummy value
else
    fseek(fid,0,-1);                                                       % set position descriptor to beginning of file
end

% read CalTable
data = fread(fid, 'double');
fclose(fid);
Nrcu = length(data) / (2 * 512);
cal = reshape(data(1:2:end) + 1i * data(2:2:end), [Nrcu, 512]);
calx = cal(1:2:end, :);
caly = cal(2:2:end, :);


