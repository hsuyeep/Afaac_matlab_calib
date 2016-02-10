% Driver code to generate images with various subband integration
% cd '/dop312_0/prasad/GPU_CORR_DAT/afaac-6/8sb_29Jan16/'
% fnames = {'SB295.vis', 'SB296.vis', 'SB297.vis', 'SB298.vis', 'SB299.vis', 'SB300.vis', 'SB301.vis', 'SB302.vis'};
cd '/home/prasad'
fnames = {'ais006dat/4100.vis', 'ais006dat/4101.vis', 'ais006dat/4102.vis', 'ais006dat/4103.vis', 'ais007dat/4104.vis', 'ais007dat/4105.vis', 'ais007dat/4106.vis', 'ais007dat/4107.vis'};

% Generate uncalibrated images from every subband over entire timerange.
obs.skip = 60;
obs.cal = 0;
for i = 1:length (fnames)
    try
        obs.sub = 294+i;
        obs.freq = obs.sub*195312.5;
        imgcorrvis ({fnames{i}}, obs, []);
    catch ME;
    end;        
end;
    
% Generate uncalibrated images: set obs.cal = 0 in imgcorrvis.m
%for i = 2:length(fnames)
%    proc_fnames = fnames (1:i);
%    imgcorrvis (proc_fnames, [], []);
%end;
