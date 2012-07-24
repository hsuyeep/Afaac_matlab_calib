function [alpha, delta] = lmtoradec(l, m, JD)

% [alpha, delta] = lmtoradec(l, m, JD)
%
% conversion from (l,m) to (ra,dec) for Dwingeloo
%
% Arguments
% l  : vector containing l-coordinates
% m  : vector containing the corresponding m-coordinates
% JD : GMT in Julian days
%
% Return values
% alpha : length(l)-by-length(m) matrix with values for the RA coordinate
%         expressed in radians. Entries corresponding to (l, m)-positions
%         that do not physically exist (sqrt(l^2 + m^2) > 1) are replaced
%         by NaN.
% delta : length(l)-by-length(m) matrix with values for the DEC coordinate
%         expressed in radians. Entries corresponding to (l, m)-positions
%         that do not physically exist (sqrt(l^2 + m^2) > 1) are replaced
%         by NaN.
%
% SJW, 2004

L = (6 + 23/60 + 41/3600);  % geographic latitude Dwingeloo in degrees
B = (52 + 48/60 + 42/3600); % geographic longitude Dwingeloo in degrees

a1 = 24110.54841;
a2 = 8640184.812866;
a3 = 0.093104;
a4 = -6.2e-6;

% Time in Julian centuries
TU = (JD - 2451545)/36525;

% Greenwich Star Time in seconds
GST = (JD+0.5)*86400 + a1 + a2*TU + a3*TU^2 + a4*TU^3;

% Local Siderderial Time in radians
% 240: number of seconds per degree
LST =  ((GST + L*240) / 240) * pi / 180;

alpha0 = LST;              % Zenith projection, in radians
delta0 = B * pi / 180;     % ibid.

% (l, m) to (alpha, delta)
for lidx = 1:length(l)
  for midx = 1:length(m)
    dist = l(lidx)^2 + m(midx)^2;
    if (dist < 1)
      alpha(lidx, midx) = atan2(l(lidx), cos(delta0) * sqrt(1 - dist) - m(midx) * sin(delta0));;
      delta(lidx, midx) = asin(sqrt(1 - dist) * sin(delta0) + m(midx) * cos(delta0));
    else
      delta(lidx, midx) = NaN;
      alpha(lidx, midx) = NaN;
    end
  end
end

% alpha on interval [0, 2 * pi>
alpha = mod(alpha0 + alpha, 2 * pi);

