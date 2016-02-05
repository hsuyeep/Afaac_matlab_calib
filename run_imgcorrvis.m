% Driver code to generate images with various subband integration
cd '/dop312_0/prasad/GPU_CORR_DAT/afaac-6/8sb_29Jan16/'
fnames = {'SB295.vis', 'SB296.vis', 'SB297.vis', 'SB298.vis', 'SB299.vis', 'SB300.vis', 'SB301.vis', 'SB302.vis'};

% Generate uncalibrated images: set obs.cal = 0 in imgcorrvis.m
for i = 2:length(fnames)
    proc_fnames = fnames (1:i);
    imgcorrvis (proc_fnames, [], []);
end;
