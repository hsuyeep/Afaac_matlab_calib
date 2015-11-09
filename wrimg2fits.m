% Script to write out AARTFAAC image data as FITS files, which can be 
% displayed in standard astronomical software like ds9.
% pep/06Nov15
% Arguments:
%   fid  :  File id of a writeable file.
%   img  :  img structure as returned by fft_imager   
%   fhdu :  FITS HDU can be presupplied for efficiency.

function wrimg2fits (fid, img, fhdu)
    
    % Determine the RA/DEC of this timeinstant based on observing time.
    [ra, dec] = convmjdsec2radec (img.tobs); % Result in rad.

    if (isempty (fhdu))
        % Generate a FITS header
        fprintf (2, '### Currently unimplemented!');
        return;
    end;
    
    

