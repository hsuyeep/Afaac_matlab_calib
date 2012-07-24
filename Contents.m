% LOFAR Station Calibration Package
%
% This software package contains the Matlab functions required to reduce LOFAR
% station calibration observations. A general overview of the station
% calibration pipeline can be found in [1], while [2] provides a step-by-step
% description of the complete station calibration process from doing the
% observation to the production of the calibration table.
%
% References
% [1] Stefan J. Wijnholds, "Overview of the Initial Production Version of the
%     Station Calibration Pipeline", Internal report LOFAR-ASTRON-MEM-263,
%     ASTRON, Dwingeloo (The Netherlands), January 2011
% [2] Stefan J. Wijnholds, "Station Calibration Cookbook", Internal report
%     LOFAR-ASTRON-MAN-049, ASTRON, Dwingeloo (The Netherlands), January 2011
%
% SJW, December 2010
%
% Calibration functions and scripts
%   calconfig          - script to set up automatic reduction of observation
%   calibrateData      - performs automatic reduction of observation
%   statcal            - wrapper to provide convenient interface to cal_ext
%   cal_ext            - performs multisource array calibration
%   findAoverT         - provides indication of A/T towards Cas A
%   makeCalTable       - script to produce calibration table
%   readCalTable       - read calibration table from calibration table file
%   writeCalTable      - saves calibration table in appropriate format
%   XYPhaseDelta       - change phase offset between x- and y-arrays
%   inspectData        - inspect subband statistics of an observation
%   acm2skyimage       - DFT imager
%   validate           - wrapper for validate_acc
%   validate_acc       - validate calibration by all-sky image comparison
%
% Support functions
%   parseAntennaArrays - extracts configuration from AntennaArrays.conf file
%   parseAntennaField  - extracts configuration from AntennaField.conf file
%   RFIdetection       - finds RFI free subbands
%   detectBadElem      - finds broken elements by comparison of pass bands
%   JulianDay          - converts UTC time to Julian Day number
%   SunRaDec           - gives position of the Sun in (ra, dec)
%   radectolm          - converts (ra, dec) coords to (l, m) coords
%   khatrirao          - computes Khatri-Rao product of two matrices
%   computeAlphaW      - estimates gain quotient matrix
%   checkCal           - performs solution based flagging
%   fitAmplModel       - script to fit amplitude model to gain solutions
%   fitPhaseModel      - script to fit phase model to gain solutions

