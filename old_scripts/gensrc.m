% Script to determine estimated source gains relative to each other, and estimated source positions wrt. 3C catalog. Operate on estimated calibration solutions using aartfaac_demo_25Jan12.
% pep/25Jan12

dirname = '../cookdat/per_snapshot_cal_sol/';
% [~, lsout] = dos(['ls -1 ' dirname '*clean.txt']);
[~, lsout] = dos(['ls -1 ' dirname '*clean_var.mat']);

datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
disp (['Found ' num2str(nfiles) ' files']);
Nelem = 288; % NOTE Hardcoded

% Needed do be done only once: load ITRF positions of antennas
fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
posITRF = zeros(Nelem, 3);
posITRF(:, 1) = posfiledata{3};
posITRF(:, 2) = posfiledata{4};
posITRF(:, 3) = posfiledata{5};

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF * rotmat;

srcsel =  [324, 283, 88, 179, 0];

% peek at the data once for matrix ranges
ccm_fname = datafiles {1}{1};
load (ccm_fname);

% Generate matrices to hold solutions across all timeslices
nsrc = size (sigmas1);
cal_tobs = zeros (nfiles);
cal1_sigmas1 = zeros (nfiles, nsrc(2));
%cal1_Sigman1 = zeros (nfiles, size (Sigman1));
cal1_lsrch  = zeros (nfiles, nsrc(2));
cal1_msrch  = zeros (nfiles, nsrc(2));

% hold off
for fnum = 1:nfiles
    disp (['Processing file ' num2str(fnum) ' of ' num2str(nfiles)]);
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);
    cal1_sigmas1 (fnum, :) = sigmas1; 
    %cal1_Sigman1 (fnum, :) = Sigman1;
    cal1_lsrch (fnum, :) = lsrchat;
    cal1_msrch (fnum, :) = msrchat;
    
    cal_tobs (fnum) = tobs;
    
   
    % figure (1);
    % plot (angle (cal1))
    % xlabel ('Dipole number');
    % ylabel ('Phase (rad)');
    % title (['Estimated per element gain phases:' datestr(tobs)]);
    
    %figure (2);
    % plot (abs (cal1));
    % xlabel ('Dipole number');
    % ylabel ('Gain');
    % title (['Estimated per element gain amplitudes:' datestr(tobs)]);
   
    %skymapcal1 = acm2skyimage(acccal1, poslocal(:, 1), poslocal(:, 2), freq, l, m);
end

x = 1:nfiles; y = 1:Nelem;
ylabel ('Estimated relative source gain');
xlabel ('Time');
plot (cal1_sigmas1 (:, 1));
