function writeCalTable(calx, caly, outfile);

% writeCalTable(calx, caly, outfile)
%
% Write calibration data to a properly formatted calibration table file.
%
% Arguments
%   calx    : Nelem x Nsb (48 x 512 for Dutch stations, 96 x 512 for European
%             stations) matrix containing a complex valued correction factor
%             for each element and each subband for the array of x-elements.
%   caly    : corresponding result for the array of y-elements
%   outfile : name of the file containing the calibration table
%
% SJW, December 2010

% initialize and fill the station calibration table
[Nant, Nfreq] = size(calx);
caltable = zeros(4 * Nant, Nfreq);
caltable(1:4:end, :) = real(calx);
caltable(2:4:end, :) = imag(calx);
caltable(3:4:end, :) = real(caly);
caltable(4:4:end, :) = imag(caly);

% write the result to the specified file
fid = fopen(outfile, 'w');
N = fwrite(fid, caltable, 'double');
if (N ~= length(caltable(:)))
    disp('The data were note properly stored.')
end
fclose(fid);
