function [clean, teststat] = RFIdetection(acc, N, RFItol)

% [clean, teststat] = RFIdeteection(acc, N, RFItol)
%
% Find RFI free subbands using the detector described in [1]
%
% arguments
% acc    : Nelem x Nelem x Nch array covariance matrix
% N      : number of samples over which the data has been integrated
% RFItol : set tolerance to this times the standard deviation
%
% return values
% clean    : Nelem x 1 selection matrix with 1 in the entries corresponding to
%            RFI free subbands and 0 otherwise
% teststat : Nelem x 1 vector containing the value of the teststatistic
%            based on which the decisions were made
%
% Reference
% [1] Stefan J. Wijnholds, "A spatial RFI detection method for LOFAR
% stations", LOFAR-ASTRON-MEM-181, rev. 1.1, July 4, 2006
%
% SJW, 2009
% modified by SJW, December 2010

nelem = size(acc, 1);
nch = size(acc, 3);
teststat = zeros(nch, 1);
treshold = RFItol * sqrt(2 * nelem / N);

for ch = 1:nch
    Rhat = squeeze(acc(:, :, ch));
    Rwhite = Rhat ./ sqrt(diag(Rhat) * diag(Rhat).');
    teststat(ch) = norm(Rwhite, 'fro').^2;
end

clean = max([abs(teststat - [teststat(2:end); 0]),  abs(teststat - [0; teststat(1:end-1)])], [], 2) < 10 * treshold;
