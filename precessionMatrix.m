function precmat = precessionMatrix(JD)

% precmat = precessionMatrix(JD)
%
% Compute rotation matrix describing the precession from J2000 to the
% specified target time in Julian days, computation according to the
% Astronomical Almanac 2008, pp B28-B31.
%
% Assumed Cartesian coordinate definitions:
% - positive x-axis of through (alpha, delta) = (0, 0)
% - positive y-axs through (alpha, delta) = (90, 0) degrees
% - positive z-axis through the NCP.
%
% argument
% JD : target time in Julian dayss, scalar
%
% return value
% precmat : 3 x 3 rotation matrix describing precession
%
% SJW & MAB, 17 May 2011

% precession time in centuries
precT = (JD - 2451545.0) / 36525.0;

% polynomial to compute coefficient psi_a in radians
polcoeff = [-0.001147, -1.07259, 5038.47875, 0].';
psi_a = (polyval(polcoeff, precT) / 3600) * pi / 180; % verified

% epsilon0 in radians
epsilon0 = (23 + 26 / 60 + 21.448 / 3600) * pi / 180; % verified

% polynomial to compute coefficient omega_a in radians
polcoeff = [-0.007726, 0.05127, -0.02524, (epsilon0 * 180 / pi) * 3600].';
omega_a = (polyval(polcoeff, precT) / 3600) * pi / 180; % verified

% polynomial to compute coefficient chi_a in radians
polcoeff = [-0.001125, -2.38064, 10.5526, 0].';
chi_a = (polyval(polcoeff, precT) / 3600) * pi / 180; % verified

% precession matrix
s1 = sin(epsilon0);
s2 = sin(-psi_a);
s3 = -sin(omega_a);
s4 = sin(chi_a);
c1 = cos(epsilon0);
c2 = cos(-psi_a);
c3 = cos(-omega_a);
c4 = cos(chi_a);
precmat = [c4 * c2 - s2 * s4 * c3, c4 * s2 * c1 + s4 * c3 * c2 * c1 - s1 * s4 * s3, c4 * s2 * s1 + s4 * c3 * c2 * s1 + c1 * s4 * s3; ...
           -s4 * c2 - s2 * c4 * c3, c4 * c3 * c2 * c1 - s4 * s2 * c1 - s1 * c4 * s3, c4 * c3 * c2 * s1 + c1 * c4 * s3 - s4 * s2 * s1; ...
           s2 * s3, -s3 * c2 * c1 - s1 * c3, c3 * c1 - s3 * c2 * s1]; % verified
