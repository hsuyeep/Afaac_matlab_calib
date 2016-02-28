% Driver code to generate images with various subband integration

% 10sec stretch from 29Jan16 data, all subbands integrated
 cd '/dop312_0/prasad/GPU_CORR_DAT/afaac-6/8sb_29Jan16/'
fnames = {'SB295.vis', 'SB296.vis', 'SB297.vis', 'SB298.vis', 'SB299.vis', 'SB300.vis', 'SB301.vis', 'SB302.vis'};
obs.sub = [295:302];
obs.flagant_x = [18, 84, 142, 168, 262];
obs.flagant_y = [18, 84, 142, 168, 262];
obs.cal = 1;
obs.stokes = 4;
%obs.imgspectint = 1;
%subs = [295:302];
%imgcorrvis (fnames, obs, []);

% SB002 run 01 on newly generated .vis visibilities from MS.
%cd /home/prasad/struisdat/fhuizing
%fnames = {'S293_C63_M9_T20131120-101100.vis'};
%obs.flagant_x = [129, 140, 149];
%obs.flagant_y = [129, 140, 149];
%obs.cal = 1;
%obs.sub = 293;
%obs.stokes = 4;
%obs.imgspectint = 1;
%obs.skip = 30;
%imgcorrvis (fnames(1), obs, '/home/prasad/struisdat/peeyush/S293_C63_M9_T20131120-101100');

% Full 10min stretch of 8 subband data, generating an image every minute.
%cd '/home/prasad'
%fnames = {'ais006dat/4100.vis', 'ais006dat/4101.vis', 'ais006dat/4102.vis', 'ais006dat/4103.vis', 'ais007dat/4104.vis', 'ais007dat/4105.vis', 'ais007dat/4106.vis', 'ais007dat/4107.vis'};
% obs.skip = 60;
%obs.flagant_x = [18, 84, 142, 168, 262];
%obs.flagant_y = [18, 84, 142, 168, 262];
%obs.cal = 1;
%subs = [295:302];
%obs.freqflag = 1;


% 10hr recording on  20Dec15
%cd '/home/prasad'
% fnames = {'ais001dat/20Dec15/SB000_1450594188.vis', 'ais003dat/20Dec15/SB002_1450594187.vis', 'ais004dat/20Dec15/SB004_1450594186.vis', 'ais005dat/20Dec15/SB006_1450594185.vis'};
%fnames = {'ais003dat/20Dec15/SB002_1450594187.vis'};
%obs.skip = 60;
%obs.sub = 297;
%obs.flagant_x = [18, 84, 142, 168, 262];
%obs.flagant_y = [18, 84, 142, 168, 262];
%obs.cal = 1;
%obs.freqflag = 0;
%obs.stokes = 4;
%obs.imgspectint = 1;
%imgcorrvis (fnames(1), obs, '/home/prasad/struisdat/peeyush/20Dec15');

% Image each subband separately
%for i = 1:length (fnames)
%    try
%        obs.sub = subs(i);
%        obs.freq = obs.sub*195312.5;
%        imgcorrvis ({fnames{i}}, obs, []);
%    catch ME;
%    end;        
%end;
    
% Combine all subbands into a single image by image plane averaging
% obs.sub = subs;
% imgcorrvis(fnames, obs, []);

% Generate a set of images with increasing spectral integration with splatting.
obs.imgspectint = 1;
obs.stokes = 4;
subs = [295:302]
for i = 1:8
    obs.sub = subs(1:i);
    imgcorrvis (fnames(1:i), obs, []);
end;

% Generate a set of images with increasing spectral integration in the image
% domain.
%obs.imgspectint = 1;
%obs.stokes = 4;
%for i = 1:8
%    obs.sub = subs(1:i);
%    imgcorrvis (fnames(1:i), obs, []);
%end;

