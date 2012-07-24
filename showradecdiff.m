% Script to show the image-plane difference between the sky projected onto RA/dec plane.
% Works on the output of genradecskyimage.m
% pep/20Apr12

load ('../cookdat/sky_radec_images_roysoc_run1.mat');

nimages = length (sky_tobs);
ragrid = linspace (0,2*pi, 1024);
decgrid = linspace (-pi/2,pi/2, 1024);

for img = 2:nimages
  % imagesc (sky_alpha*180/pi, sky_delta*180/pi, sky_radecimage_subAteam (:, :, img));
  imagesc (ragrid*180/pi, decgrid*180/pi, sky_radecimage_subAteam (:, :, img));
  colorbar;
  caxis ([0 6000]);
  set(gca, 'FontSize', 16);
  title (['All sky: ' datestr(sky_tobs(img,1))]);
  xlabel ('RA (deg)');
  ylabel ('Dec. (deg)');
  saveas (gcf, datestr(sky_tobs(img,1)), 'png');

  imdiff = sky_radecimage_subAteam (:,:,img) - sky_radecimage_subAteam (:,:,img-1);
  % imagesc (sky_alpha*180/pi, sky_delta*180/pi, imdiff);
  imagesc (ragrid*180/pi, decgrid*180/pi, imdiff);
  colorbar;
  caxis ([0 6000]);
  set(gca, 'FontSize', 16);
  title (['Difference image: ' datestr(sky_tobs(img,1))]);
  xlabel ('RA (deg)');
  ylabel ('Dec. (deg)');
  saveas (gcf, strcat (datestr(sky_tobs(img,1)), '_diff'), 'png');
end
