function cal = checkCal(cal, sel, Nch, nelem)

% cal = checkCal(cal, sel, Nch, nelem)
%
% Remove bad calibration solutions and produce a gain correction matrix of
% appropriate size (Nch x Nelem)
%
% arguments
% cal   : Nch x Ngoodelem matrix with raw calibration result
% sel   : Nelem x 1 vector indicating the good elements used for calibration
% Nch   : number of frequency channels
% Nelem : number of elements
%
% return value
% cal : Nch x Nelem matrix with filtered calibration data
%
% SJW, 2009

% make sure that the median calibration correction equals 1
cal = cal ./ repmat(median(abs(cal), 2), [1 size(cal, 2)]);
% make sure blanked channels also get a well defined value
cal(isnan(cal)) = 0;
% check for too large deviations between subbands
delta = max(abs(cal - [cal(2:end, :); zeros(1, size(cal, 2))]), abs(cal - [zeros(1, size(cal, 2)); cal(1:end-1, :)]));
selcalf = sum(delta < 0.3, 2) > 0.5 * nelem;
% check for deviating antennas
selcalant = sum(abs(cal(selcalf, :)) < 3) == max(sum(abs(cal(selcalf, :)) < 3));
% blank ill subbands and antennas
cal(~selcalf, :) = 0;
cal(:, ~selcalant) = 0;
caltemp = zeros(Nch, nelem);
caltemp(:, sel) = cal;
cal = caltemp;
