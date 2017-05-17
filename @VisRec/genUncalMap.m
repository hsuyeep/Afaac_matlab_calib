% Function to generate an uncalibrated map from the observed
% visibilities. Returns generated images in the img structure
% Both XX and YY images are returned.
%
% Arguments:
%   arrayconfig : One of 'lba_outer', 'lba_inner'
% Returns:
%   img : structure containing generated images.
function img = genUncalMap (obj, arrayconfig)
    assert (obj.freq ~= 0);
    load ('poslocal_outer.mat', 'poslocal');
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	gparm.type = 'pillbox';
    gparm.lim  = 0;
    gparm.duv = 0.5;				% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    gparm.Nuv = 512;				% size of gridded visibility matrix
    gparm.uvpad = 512;				% specifies if any padding needs to be added
    gparm.fft  = 1;

    img.tobs = obj.trecstart;
    img.freq = obj.freq;
    [~, img.map_xx, ~, img.l, img.m] = ... 
        fft_imager_sjw_radec (obj.acm_xx(:), uloc(:), vloc(:), ... 
			gparm, [], [], img.tobs, img.freq, 0);
    [~, img.map_yy, ~, img.l, img.m] = ... 
        fft_imager_sjw_radec (obj.acm_yy(:), uloc(:), vloc(:), ... 
			gparm, [], [], img.tobs, img.freq, 0);
end;
