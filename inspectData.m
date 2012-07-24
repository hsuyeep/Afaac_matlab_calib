function inspectData(dirname)

% inspectData(dirname)
%
% Inspect subband statistics collected during monitoring observations.
% This function composes a list of all subband statistics files in the
% specified directory and plots the data in each file consecutively using
% the RCU number as plot title.
%
% SJW, December 2010

% compose list of all subband statistics files in specified directory
[status, lsout] = dos(['ls -1 ' dirname '*sst*']);
datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});

% loop over all data files showing their contents
for idx = 1:nfiles
    fid = fopen(datafiles{1}{idx}, 'r');
    data = fread(fid, 'double');
    fclose(fid);
    tfdata = reshape(data, [512, length(data)/512]);
    imagesc(log10(tfdata));
    title(datafiles{1}{idx}(length(dirname)+21:length(dirname)+26));
    pause
end

