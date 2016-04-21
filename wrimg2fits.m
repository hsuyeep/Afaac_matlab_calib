% Script to write out AARTFAAC image data as FITS files, which can be 
% displayed in standard astronomical software like ds9.
% pep/06Nov15
% Arguments:
%   img  :  img structure as returned by fft_imager   
%   fname:  Name with which the file should be saved.

function wrimg2fits (img, fname);
    
    assert (all (isfield (img, {'pix2laxis', 'pix2maxis', 'tobs', 'freq', 'df', 'map'})));
    import matlab.io.*; % For low level access to the FITS routines.
    % Determine the RA/DEC of this timeinstant based on observing time.
    [ra, dec] = convmjdsectoradec (img.tobs); % Result in rad.

    try
        fid = fits.createFile (fname);
    catch ME
        ME1 = MException ('wrimg2fits:filecreationerror', 'createFile failed!');
        ME1 = addCause (ME1, ME);
        throw (ME1);
    end;
    fits.createImg (fid, 'float_img', [img.pix2laxis, img.pix2maxis, 1, 1]); 
    fits.writeKey  (fid, 'BSCALE', 1);
    fits.writeKey  (fid, 'BZERO',  0);        
    fits.writeKey  (fid, 'BMAJ' ,  1);         
    fits.writeKey  (fid, 'BMIN' ,  1);         
    fits.writeKey  (fid, 'BPA'  ,  0);         
    fits.writeKey  (fid, 'BTYPE','Intensity');
    fits.writeKey  (fid, 'OBJECT','Aartfaac image'); 
    fits.writeKey  (fid, 'BUNIT'  ,'Jy/beam');                          
    fits.writeKey  (fid, 'EQUINOX',                       2000);      
    fits.writeKey  (fid, 'RADESYS',                      'FK5');                             
    fits.writeKey  (fid, 'LONPOLE',                        180);      
    fits.writeKey  (fid, 'LATPOLE',                     52.915);      
    fits.writeKey  (fid, 'PC01_01',                          1);      
    fits.writeKey  (fid, 'PC02_01',                          0);      
    fits.writeKey  (fid, 'PC03_01',                          0);      
    fits.writeKey  (fid, 'PC04_01',                          0);      
    fits.writeKey  (fid, 'PC01_02',                          0);      
    fits.writeKey  (fid, 'PC02_02',                          1);      
    fits.writeKey  (fid, 'PC03_02',                          0);      
    fits.writeKey  (fid, 'PC04_02',                          0);      
    fits.writeKey  (fid, 'PC01_03',                          0);      
    fits.writeKey  (fid, 'PC02_03',                          0);      
    fits.writeKey  (fid, 'PC03_03',                          1);      
    fits.writeKey  (fid, 'PC04_03',                          0);      
    fits.writeKey  (fid, 'PC01_04',                          0);      
    fits.writeKey  (fid, 'PC02_04',                          0);      
    fits.writeKey  (fid, 'PC03_04',                          0);      
    fits.writeKey  (fid, 'PC04_04',                          1);      
    fits.writeKey  (fid, 'CTYPE1' ,                 'RA---SIN');                        
    fits.writeKey  (fid, 'CRVAL1' ,                  ra*180/pi);      
    fits.writeKey  (fid, 'CDELT1' ,-asind(1/(img.pix2laxis/2)));      
    fits.writeKey  (fid, 'CRPIX1' ,          img.pix2laxis/2+1);      
    fits.writeKey  (fid, 'CUNIT1' ,                      'deg');                             
    fits.writeKey  (fid, 'CTYPE2' ,                 'DEC--SIN');                        
    fits.writeKey  (fid, 'CRVAL2' ,                     52.915);      
    fits.writeKey  (fid, 'CDELT2' , asind(1/(img.pix2laxis/2)));
    fits.writeKey  (fid, 'CRPIX2' ,          img.pix2maxis/2+1);      
    fits.writeKey  (fid, 'CUNIT2' ,                      'deg');                     
    fits.writeKey  (fid, 'CTYPE3' ,                     'FREQ');                            
    fits.writeKey  (fid, 'CRVAL3' ,                   img.freq);     
    fits.writeKey  (fid, 'CDELT3' ,                     img.df);     
    fits.writeKey  (fid, 'CRPIX3' ,                         1 );     
    fits.writeKey  (fid, 'CUNIT3' ,                      'Hz' );                             
    fits.writeKey  (fid, 'CTYPE4' ,                  'STOKES' );                         
    fits.writeKey  (fid, 'CRVAL4' ,                         1 );     
    fits.writeKey  (fid, 'CDELT4' ,                         1 );     
    fits.writeKey  (fid, 'CRPIX4' ,                         1 );     
    fits.writeKey  (fid, 'CUNIT4' ,    'stokes-unit'          );           
    fits.writeKey  (fid, 'PV2_1'  ,                         0 );     
    fits.writeKey  (fid, 'PV2_2'  ,                         0 );     
    fits.writeKey  (fid, 'RESTFRQ',                   img.freq);
    fits.writeKey  (fid, 'SPECSYS',                    'LSRK' );                   
    fits.writeKey  (fid, 'ALTRVAL',                         0);
    fits.writeKey  (fid, 'ALTRPIX',                         1);
    fits.writeKey  (fid, 'VELREF' ,                       257);
    fits.writeKey  (fid, 'TELESCOP',    'AARTFAAC'            )             
    fits.writeKey  (fid, 'OBSERVER',   'AARTFAAC Project'     )            
    fits.writeKey  (fid, 'DATE-OBS',datestr(mjdsec2datenum(img.tobs), 'yyyy-mm-ddTHH:MM:SS.FFF'));
    fits.writeKey  (fid, 'TIMESYS',     'UTC'                 );              
    fits.writeKey  (fid, 'OBSRA'  ,                     ra*180/pi);      
    fits.writeKey  (fid, 'OBSDEC' ,                     52.915);      
    fits.writeKey  (fid, 'OBSGEO-X',                3.8266e+06);% CS002 center ITRF location    
    fits.writeKey  (fid, 'OBSGEO-Y',                4.6102e+05);      
    fits.writeKey  (fid, 'OBSGEO-Z',                5.0649e+06);      
    fits.writeKey  (fid, 'DATE'   , datestr(now, 'yyyy-mm-ddTHH:MM:SS.FFF'));
    fits.writeKey  (fid, 'ORIGIN' ,  'wrimg2fits.m');

    % fliplr added by pep/14Apr16 to ensure proper coordinates on the generated
    % images.
    fits.writeImg  (fid, fliplr(img.map)); 
    fits.closeFile (fid);
