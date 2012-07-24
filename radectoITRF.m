function ITRFpos = radectoITRF(ra, dec, epoch, JD)

% ITRFpos = radectoITRF(ra, dec, epoch, JD)
%
% Conversion from source positions in right ascension and declination to a
% Cartesian direction vector in ITRF coordinates.
%
% arguments
% ra    : Nsrc x 1 vector with right ascension for each source in rad in
%         either B1950 or J2000 coordinates
% dec   : Nsrc x 1 vector with declination for each source in rad in either
%         B1950 or J2000 coordinates
% epoch : Nsrc x 1 vector with boolean value, true if the corresponding
%         source position is given in B1950 coordinates, false for J2000
% JD    : target time in Julian days
%
% return value
% ITRFpos : Nsrc x 3 matrix with direction vector of unit lenght for each
%           source
%
% SJW & MAB, 17 May 2011

% convert (ra, dec) to Cartesian coordinates
Nsrc = length(ra);
radecCart = zeros(Nsrc, 3);
[radecCart(:, 1), radecCart(:, 2), radecCart(:, 3)] = sph2cart(ra, dec, 1);

% convert B1950 to J2000
precMat = precessionMatrix(JulianDay(datenum(1950, 1, 1, 0, 0, 0)));
radecCart(epoch, :) = radecCart(epoch, :) * precMat; % verified

% compute apparent positions
precMat = precessionMatrix(JD);
app_pos = radecCart * precMat.'; % verified

% Greenwich mean sidereal time in seconds
a1 = 24110.54841;
a2 = 8640184.812866;
a3 = 0.093104;
a4 = -6.2e-6;
polcoeff = [a4, a3, a2, a1];
TU = (floor(JD) + 0.5 - 2451545) / 36525;
GMST = (JD - floor(JD) - 0.5) * 86400 * 1.002737811906 + polyval(polcoeff, TU);
% conversion to radians
GMST = (GMST / 86400) * 2 * pi; % verified

% Convert apparent positions to ITRF positions
rotmat = [ cos(GMST), sin(GMST), 0; ...
          -sin(GMST), cos(GMST), 0; ...
                   0,         0, 1];
ITRFpos = app_pos * rotmat.';
