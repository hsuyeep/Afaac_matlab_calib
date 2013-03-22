function validate(datafile, calfile, posfile, rcumode, sbsel)

% validate(datafile, calfile, posfile, rcumode, sbsel)
%
% Validation of the calibration table using data from the calibration
% observation by comparison of images before and after calibration. The
% images produced are averaged over the use specified subband selection.
% This function is a wrapper function of validate_acc.
%
% Arguemens
% datafile : acc-file from calibration observation
% calfile  : calibration table file as produced by writeCalTable
% posfile  : AntennaField.conf file of appropriate station
% rcumode  : RCU mode number to validate
% sbsel    : subband selection to use for validation
%
% SJW, October 2011

% read the data
fid = fopen(datafile, 'r');
data = fread(fid, 'double');
fclose(fid);
Nsb = 512;
Nrcu = sqrt(length(data) / (2 * Nsb));
acc = reshape(data(1:2:end) + 1i * data(2:2:end), [Nrcu, Nrcu, Nsb]);

% call validate_acc to make the images
refidx = strfind(datafile, '_acc_');
tobs = datenum(datafile(refidx-15:refidx-1), 'yyyymmdd_HHMMss');
validate_acc(acc(:, :, sbsel), calfile, posfile, rcumode, sbsel, tobs);
