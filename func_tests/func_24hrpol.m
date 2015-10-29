% Script to generate stokes I and stokes Q images from XX and YY 24hr. images.
% The output is written out as fits files without a WCS.
% Note that since regenerating the image data is timeconsuming, we operate
% directly on the resulting .fig files.
% pep/11May15

clear all; close  all;
% Load calibrated and A-team subtracted data (already generated).
load ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_24hrobs/24hr_radec_cal_ateamsub_1secaccum.fig', '-mat');
XX_imgdat = hgS_070000.children(1).children(1).properties.CData;
clear hgS_070000;

load ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_24hrobs/24hr_FFTimg_radec_cal_ateamsub_1secper1min_accum_YY.fig', '-mat');
YY_imgdat = hgS_070000.children(1).children(1).properties.CData;

% Generate Stokes I
I_imgdat = XX_imgdat + YY_imgdat;
Q_imgdat = XX_imgdat - YY_imgdat;

% Display all images
al_lin = linspace (0, 2*pi, 512);
de_lin = linspace (-pi/2, pi/2, 512);

figure
imagesc (al_lin*12/pi, de_lin*180/pi, XX_imgdat); 
title ('All-sky XX image, LBA\_OUTER, 20Nov2013, 57 MHz');
xlabel ('RA (Hrs)')
ylabel ('Dec (deg)');
colorbar;
print ('all_sky_XX.eps', '-depsc');


figure
imagesc (al_lin*12/pi, de_lin*180/pi, YY_imgdat); 
title ('All-sky YY image, LBA\_OUTER, 20Nov2013, 57 MHz');
xlabel ('RA (Hrs)')
ylabel ('Dec (deg)');
colorbar;
print ('all_sky_YY.eps', '-depsc');

figure
imagesc (al_lin*12/pi, de_lin*180/pi, I_imgdat); 
title ('All-sky Stokes-I image, LBA\_OUTER, 20Nov2013, 57 MHz');
xlabel ('RA (Hrs)')
ylabel ('Dec (deg)');
colorbar;
print ('all_sky_I.eps', '-depsc');

figure
imagesc (al_lin*12/pi, de_lin*180/pi, Q_imgdat); 
title ('All-sky Stokes-Q image, LBA\_OUTER, 20Nov2013, 57 MHz');
xlabel ('RA (Hrs)')
ylabel ('Dec (deg)');
colorbar;
print ('all_sky_Q.eps', '-depsc');


% Write out as FITS.
fitswrite (XX_imgdat, 'allsky_ateamsub_XX.fits');
fitswrite (YY_imgdat, 'allsky_YY.fits');
fitswrite (I_imgdat, 'allsky_I.fits');
fitswrite (Q_imgdat, 'allsky_Q.fits');

