% calconfig
%
% This script sets up the automatic reduction of a station calibration
% observation by filling in the config struct required by the calibrateData
% function.
%
% For convenience, a string "cmd" is produced, such that the station
% calibration pipeline can be started with "eval(cmd)".
%
% SJW, December 2010

% ensure a clean workspace
clear all
close all

% directory containing the calibration data, ensure that it ends with "/"
config.dirname = '/dop147_0/wijnholds/data/caldata/CS002_rcumode_1_20110730/';

% start and stop frequency (in Hz) of band pass filter
config.bandstart = 10e6;
config.bandstop = 90e6;

% start and stop frequency (in Hz) of Nyqyuist zone
config.startfreq = 0e6;
config.stopfreq = 100e6;

% source selection (index in 3C catalog), use index 0 for the Sun
config.srcsel = [324, 283, 88, 179, 0];

% AntennaArrays.conf file and name of array configuration to be used
config.posfile = 'AntennaFields/AntennaFieldCS002.conf';
config.rcumode = 1;
config.EU = 0;

% baseline specific flags to filter, e.g., crosstalk
config.uvmask = eye(48); %kron(eye(3), ones(16));

% set RFI tolereance in standard deviations
% typical values: 1 for LBA, 5 for HBA
config.RFItol = 1;

% name of the file in which the results will be stored
config.outfile = 'CS002_test20110928.mat';

% the command string cmd
cmd = 'calibrateData(config)';
