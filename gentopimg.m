% Script to generate an image for all channels and times in a subband from
% the TOP project CygA observation.
% pep/13Dec12

fid = fopen ('~/WORK/TOP_CygA/L77781_SAP000_SB123_uv.bin', 'rb');
nchan = 3; nant = 34; % Hardcoded
load 'posITRF_L77781.mat';

% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF_L77781 * rotmat;

hdl = zeros (1, 3);
hdl(1) = figure; hdl(2) = figure; hdl(3) = figure;
l = -1:0.01:1; m = l;
for ind = 1:100
	for ch=1:3
		[acc, tobs, freq] = readms2float (fid, -1, -1, nant);
		skymap = acm2skyimage (acc, posITRF_L77781(:,1), posITRF_L77781(:,2), freq, l, m);
		figure (hdl(ch));	
		imagesc (skymap); 
		colorbar;
	end;
	pause;
end;
	
