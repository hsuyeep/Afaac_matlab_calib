function skymap = acm2skyimage(acm, xpos, ypos, zpos, freq, l, m)

% skymap = acm2skyimag(acm, xpos, ypos, freq, l, m)
%
% conversion of ACM to sky image
%
% arguments
% acm    : nelem x nelem x nchannel array correlation matrix
% xpos   : x-position of antennas (vector of length nelem)
% ypos   : y-position of antennas (vector of length nelem)
% freq   : frequencies in Hz (vector of length nchannel)
% l, m   : l and m points to which to direct the array (vectors)
%
% return value
% skymap : length(l) x length(m) x nchannel matrix containing the
%          resulting sky maps
%
% SJW, 2004
%
% modified April 19, 2005 by SJW: correct coordinate conventions
% modified July 20, 2006 by SJW: optimization for release in LOFAR package

[nelem, dummy, nchannel] = size(acm);

c = 2.9979245e8;

skymap = zeros(length(l), length(m), nchannel);

for nch = 1:nchannel
  disp(['acm2skyimage: working on channel ' num2str(nch) ' of ' num2str(nchannel)]);
  lambda = c / freq(nch);
  k = 2 * pi / lambda;
  wx = exp(-i * k * xpos(:) * l(:).');
  wy = exp(-i * k * ypos(:) * m(:).');
  wz = exp(-i * k * zpos(:) * (sqrt (1 - l(:).^2 - m(:).^2)).');
  for lidx = 1:length(l)
    for midx = 1:length(m)
      weight = wx(:, lidx) .* wy(:, midx) .* wz(:,midx);
      skymap(lidx, midx, nch) = real(weight' * (acm(:, :, nch)) * weight);
    end
  end
end
