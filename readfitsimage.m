% Script to read FITS images, for interfacing with the AARTFAAC calim pipeline.
% Note: The assumption is that there is only one image per file, as generated
% by conversion of casacore images to FITS.
% pep/04Nov15

% Arguments:
%   fname:  FITS File name.
%   deb  :  Bool to generate verbose messages. 

% Returns:
%	img.tobs     :   The MJD sec. time corresponding to this map.
%	img.freq     :   The Freq. in Hz. corresponding to this map.
% 	img.pix2laxis:   Number of pixels along the l-axis.
% 	img.pix2maxis:   Number of pixels along the m-axis.
%	img.map      :   Matrix containing the (square) image of size (pix2laxis X pix2maxis)
%			         to be written out. This is always real.
%   img.l,m      :   Vector containing the l and m coordinates corresponding to 
%				     generated image.


function [img] = readfitsimage (fname, deb)
    assert (isempty (fname) == 0)
    info = fitsinfo (fname); 

    ind = find (strcmp (info.PrimaryData.Keywords, 'NAXIS1'));
    img.pix2laxis = info.PrimaryData.Keywords{ind,2}; % Assuming NAXIS1 is l, validate!

    ind = find (strcmp (info.PrimaryData.Keywords, 'NAXIS2'));
    img.pix2maxis = info.PrimaryData.Keywords{ind,2}; % Assuming NAXIS2 is m, validate!

    ind = find (strcmp (info.PrimaryData.Keywords, 'CRVAL3'));
    img.freq = info.PrimaryData.Keywords{ind,2}; % Assuming CTYPE3 is freq, validate!

    ind = find (strcmp (info.PrimaryData.Keywords, 'DATE-OBS'));
    img.tobs = datenum2mjdsec (datenum (info.PrimaryData.Keywords{ind,2}, 'yyyy-mm-ddTHH:MM:SS'));
    
    ind = find (strcmp (info.PrimaryData.Keywords, 'CDELT1'));
    extent = sin (info.PrimaryData.Keywords{ind,2} * img.pix2laxis);
    img.l = linspace (-extent, extent, img.pix2laxis);

    ind = find (strcmp (info.PrimaryData.Keywords, 'CDELT2'));
    extent = sin (info.PrimaryData.Keywords{ind,2} * img.pix2laxis);
    img.m = linspace (-extent, extent, img.pix2maxis);

    img.map = fitsread (fname);
