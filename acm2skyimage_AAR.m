function [skymap_num,skymap_den,skymap_classical] = acm2skyimage(Tacm, xpos, ypos, freq, l, m)

% skymap = acm2skyimag(Tacm, xpos, ypos, freq, l, m)
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
%
% Modified by Amir Leshem 2012


[nelem, dummy, nchannel] = size(Tacm);

c = 2.9979245e8;

skymap_classical = zeros(length(l), length(m),nchannel);
skymap_num = zeros(length(l), length(m),nchannel);
skymap_den=  zeros(length(l), length(m),nchannel);
for nch = 1:nchannel
  disp(['acm2skyimage: working on channel ' num2str(nch) ' of ' num2str(nchannel)]);
  lambda = c / freq(nch);
  k = 2 * pi / lambda;
  wx = exp(-i * k * xpos(:) * l(:).');
  wy = exp(-i * k * ypos(:) * m(:).'); 
  acm=Tacm(:,:,nch);
  invacm=inv(acm(:, :, nch));
  invacm2=invacm^2;
  for lidx = 1:length(l)
      if lidx==20*floor(lidx/20)
      disp(['l=',int2str(lidx),'/',int2str(length(l))]);
      end
    for midx = 1:length(m)
      weight = wx(:, lidx) .* wy(:, midx);
      skymap_classical(lidx, midx, nch) = real(weight' * acm* weight);
      skymap_num(lidx, midx, nch) = real(weight' * invacm* weight);
      skymap_den(lidx, midx, nch) = real(weight' * invacm2* weight);
      end
  end
end
clear Tacm acm invacm