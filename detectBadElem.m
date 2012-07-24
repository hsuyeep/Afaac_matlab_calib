function [sel, ac] = detectBadElem(acc)

% [sel, ac] = detectBadElem(acc)
%
% Find broken elements based on a comparison of the pass bands
%
% argument
% acc : Nelem x Nelem x Nch array covariance matrix
%
% return values
% sel : Nelem x 1 vector containing a 1 at the entries correspoding to
%       properly working elements and a 0 at the entries corresponding to
%       broken elements
% ac  : Nelem x Nch autocorrelation matrix containing the autocorrelation
%       spectra based on which the decisions were made
%
% SJW, 2009

nelem = size(acc, 1);
nch = size(acc, 3);
ac = zeros(nelem, nch);
maxdev = 2; % in dB
for idx = 1:nelem
    ac(idx, :) = squeeze(acc(idx, idx, :));
end
selfactor = ac ./ repmat(median(ac), nelem, 1);
sel = sum((selfactor < 10^(maxdev/10)).' & (1 ./ selfactor < 10^(maxdev/10)).') == nch;
