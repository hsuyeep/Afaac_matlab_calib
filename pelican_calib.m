% Script to carry out calibration to convergence of a given ACM, for a given 
% timeslice and frequency bin.
function [calvis, gainsol, sigmas, sigman, good] = pelican_calib (acc, t_obs, freq, uvflag)
%
% Arguments:
% t_obs : Time in MJD secs, as extracted from a MS
% acc   : Filled, square ACM of size 288x288 complex numbers
% freq  : The frequency of the observation in Hz.
% uvflag: Mask for flagged visibilities, which are ignored 

% Return values:
% calvis : Calibrated visibilities, ready for imaging
% gainsol: Estimated per-element complex gains
% sigmas : Estimated flux of model sources
% sigman : Estimated system temprature per baseline
% good   : Boolean encoding goodness of calibration solution
% Dependencies:
% statcal_vlss.m: Carry out the actual calibration
%  --> SunRaDec.m, radectoITRF.m, khatrirao.m, 
%  --> cal_ext.m
%      --> computeAlphaW.m

    % Load once only
    persistent first_call;
    persistent srclist3CR;
    persistent posITRF;
    persistent poslocal;

    if (isempty (first_call))
        disp ('   --- Non Tracking calibration to convergence ---   ');
        disp ('Loading 3CR catalog and local antenna positions.');
    	load srclist3CR
    	load ('poslocal.mat', 'posITRF', 'poslocal'); 
        first_call = 0;
    end
    C = 299792458;               % speed of light, m/s
    lon = 6.869837540;           % longitude of CS002 in degrees
    lat = 52.915122495;          % latitude of CS002 in degrees
    restriction = 4;             % Avoid uv values with <10 wavelengths length
    maxrestriction = 60;         % Avoid uv values with >60 wavelengths length
    Nelem = 288;                 % NOTE: Hardcoded!
    nblines = Nelem * (Nelem + 1)/2; 
    normal = [0.598753, 0.072099, 0.797682].';
    srcsel =  [324, 283, 88, 179, 0]; % A team from 3CR catalog
    good = 1;                    % Optimistic approach!
	
    if (isempty (acc) | isempty (t_obs))
        disp ('ACM or time missing!');   
        return;
    end
    t_obs  = t_obs/86400 + 2400000.5; % Convert to MJD day units
    catalog = srclist3CR;
    
    [cal1, sigmas, sigman] = statcal_vlss(acc, t_obs, freq, posITRF, srcsel, normal, 10, maxrestriction, eye(Nelem), catalog);
    calvis = (cal1' * cal1) .* acc;
    gainsol = cal1;

    % Check for goodness of calibration (TODO)
    
    % Check for sun's model fitting (TODO)
    
    % WSF source position estimation of A-team (TODO)
    
    % Removal of bright sources (TODO)
    
    % Recalibrate after bright source removal (TODO)
    
    % disp (['RA of skymodel sources: ']);
    % disp ([srclist3CR(srcsel ~= 0).alpha]);
