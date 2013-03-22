% Script to generate skyimages from calibrated visibilities.
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


% create 3D matrix to store images in time order
ccm_fname = datafiles {1}{1};
load (ccm_fname);
 
z = size (l);
sky_cal2 = zeros (z(2), z(2), nfiles);
sky_subAteam = zeros (z(2), z(2), nfiles);
sky_tobs = zeros (nfiles);
sky_l = zeros (z(2), nfiles);
sky_m = zeros (z(2), nfiles);

for fnum = 1:nfiles
    disp (['Processing file ' num2str(fnum) ' of ' num2str(nfiles)]);
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);
    
    mask = NaN(length(l));
    mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;

    sky_tobs (fnum) = tobs;
    sky_l (:, fnum)  = l;
    sky_m (:, fnum) = m;
    sky_cal2(:, :, fnum) = acm2skyimage(acccal2, poslocal(:, 1), poslocal(:, 2), freq, l, m);
    % sky_subAteam (:, :, fnum) = acm2skyimage(accsubAteam, poslocal(:, 1), poslocal(:, 2), freq, l, m);
    
    
%     imagesc(l, m, skymapcal1 .* mask);
%     set(gca, 'FontSize', 16);
%     % title('extended emission, Sun and A-team removed');
%     axis equal
%     axis tight
%     xlabel('South \leftarrow m \rightarrow North');
%     ylabel('East \leftarrow l \rightarrow West');
%     set(colorbar, 'FontSize', 16);
end
outfilename = 'sky_cal2_images.mat';
save (outfilename, 'sky_cal2', 'sky_tobs', 'sky_l', 'sky_m');

% outfilename = 'sky_subAteam_images.mat';
%save (outfilename, 'sky_subAteam', 'sky_tobs', 'sky_l', 'sky_m');
