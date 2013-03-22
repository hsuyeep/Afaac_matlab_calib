% Script to operate on estimated calibration solutions using aartfaac_demo_25Jan12.
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

% Generate matrices to hold solutions across all timeslices
cal_tobs = zeros (nfiles);
cal1_ph = zeros (nfiles, Nelem);
cal1_phdiff = zeros (nfiles - 1, Nelem);

cal2_ph = zeros (nfiles, Nelem);
cal2_phdiff = zeros (nfiles - 1, Nelem);

% hold off
for fnum = 1:nfiles
    disp (['Processing file ' num2str(fnum) ' of ' num2str(nfiles)]);
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);
    cal1_ph (fnum, :) = angle (cal1); 
    cal2_ph (fnum, :) = angle (cal2);
    cal_tobs (fnum) = tobs;
    
    % For plotting phase differences from one time instant to the next
    if fnum > 2
        cal1_phdiff (fnum-1, :) = cal1_ph (fnum, :) - cal1_ph (fnum-1, :);
        cal2_phdiff (fnum-1, :) = cal2_ph (fnum, :) - cal2_ph (fnum-1, :);
    end
  
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
xlabel ('Dipole number');
ylabel ('Time');
zlabel ('Phase (rad)');
surf (x, y, cal1_ph');
shading flat;

figure 
surf (x, y, cal2_ph');
shading flat

figure
surf (1:nfiles-1, y, cal1_phdiff');
shading flat

figure
surf (1:nfiles-1, y, cal2_phdiff');
shading flat
