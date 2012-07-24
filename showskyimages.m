% Script to display images from saved .mat files containing image outputs.
% pep/25Jan12

% load 'sky_subAteam_images.mat'
load 'sky_cal2_images.mat';
%load 'sky_cal1_images.mat';


a = size (sky_l);
disp (['Found ' num2str(a(2)) ' images']);

for im = 1:a(2)
    figure;
    
    mask = NaN(a(1));
    mask(meshgrid(sky_l (:,im)).^2 + meshgrid(sky_m(:,im)).'.^2 < 1) = 1;
    
    disp (['Showing image ' num2str(im) ' of ' num2str(a(2))]);
    
    % NOTE: For displaying from skycal1_images.mat
    % imagesc(sky_l(:,im), sky_m(:,im), sky_cal1 (:,:,im) .* mask);
    
    % NOTE: For displaying from sky_cal2_images.mat
    imagesc(sky_l(:,im), sky_m(:,im), ((sky_cal2 (:,:,im).' .* mask)));
    
    
    % NOTE: For displaying from sky_subAteam_images.mat
    % imagesc(sky_l(:,im), sky_m(:,im), sky_subAteam (:,:,im).' .* mask);
    
    
    set(gca, 'FontSize', 16);
    
    % title (['extended emission, Sun and A-team removed' datestr(sky_tobs(im))]);
    title (['Station calibration' datestr(sky_tobs(im))]);
    axis equal
    axis tight
    set (gca, 'YDir', 'Normal'); % To match orientation with station images
    set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    ylabel('South \leftarrow m \rightarrow North');
    xlabel('East \leftarrow l \rightarrow West');
    
    
    set(colorbar, 'FontSize', 16);
    caxis ([-0.5e3 1.5e4]); % for calibrated images without Ateam/sun
    % removal
    % caxis ([-500 5000]); % for calibrated images with Ateam/sun removed
    pause (1);
end
