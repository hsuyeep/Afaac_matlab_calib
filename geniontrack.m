% Script to generate plots of catalog Vs. WSF positions.
% Works on the output of wsfsrcpos.m
% pep/20Apr12

fname = '../cookdat/SB000_ch30-35_5sec_3hr_wsfthphipos.mat'

load (fname);

rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];

% Convert WSF and catalog th/phi to ITRF
srcposwsf_x = cos (thsrc_wsf).* cos(phisrc_wsf);
srcposwsf_y = sin (phisrc_wsf).* cos(thsrc_wsf);
srcposwsf_z = sin(thsrc_wsf);

srcposcat_x = cos (thsrc_cat).* cos(phisrc_cat);
srcposcat_y = sin(phisrc_cat) .* cos(thsrc_wsf);
srcposcat_z = sin(thsrc_cat);

% Collect xyz coordinates for individual sources together
srcwsf1_xyz = [srcposwsf_x(1,:)' srcposwsf_y(1,:)' srcposwsf_z(1,:)'];
srcwsf2_xyz = [srcposwsf_x(2,:)' srcposwsf_y(2,:)' srcposwsf_z(2,:)'];
srcwsf3_xyz = [srcposwsf_x(3,:)' srcposwsf_y(3,:)' srcposwsf_z(3,:)'];

srccat1_xyz = [srcposcat_x(1,:)' srcposcat_y(1,:)' srcposcat_z(1,:)'];
srccat2_xyz = [srcposcat_x(2,:)' srcposcat_y(2,:)' srcposcat_z(2,:)'];
srccat3_xyz = [srcposcat_x(3,:)' srcposcat_y(3,:)' srcposcat_z(3,:)'];

% Rotate to CS002
srcwsf1_rot = rotmat' * srcwsf1_xyz';
srcwsf2_rot = rotmat' * srcwsf2_xyz';
srcwsf3_rot = rotmat' * srcwsf3_xyz';

srccat1_rot = rotmat' * srccat1_xyz';
srccat2_rot = rotmat' * srccat2_xyz';
srccat3_rot = rotmat' * srccat3_xyz';

% Calculate ele. and azi.
srccat1_el = asin (srccat1_rot(3,:));
srccat1_azi = atan2 (srccat1_rot(1,:), srccat1_rot (2,:));
srccat2_el = asin (srccat2_rot(3,:));
srccat2_azi = atan2 (srccat2_rot(1,:), srccat2_rot (2,:));
srccat3_el = asin (srccat3_rot(3,:));
srccat3_azi = atan2 (srccat3_rot(1,:), srccat3_rot (2,:));

srcwsf1_el = asin (srcwsf1_rot(3,:));
srcwsf1_azi = atan2 (srcwsf1_rot(1,:), srcwsf1_rot (2,:));
srcwsf2_el = asin (srcwsf2_rot(3,:));
srcwsf2_azi = atan2 (srcwsf2_rot(1,:), srcwsf2_rot (2,:));
srcwsf3_el = asin (srcwsf3_rot(3,:));
srcwsf3_azi = atan2 (srcwsf3_rot(1,:), srcwsf3_rot (2,:));

% Plot the tracks, note that there is a -pi to pi transition in srcwsf3 (1)
%figure;
%plot (srccat1_azi*180/pi, srccat1_el*180/pi);
%hold on
%plot (srcwsf1_azi*180/pi, srcwsf1_el*180/pi, 'r');
%set(gca, 'FontSize', 16);
%title ('Predicted Vs. observed track of CasA over 3 Hrs');
%xlabel ('Azimuth (deg)');
%ylabel ('Elevation (deg)');
%
%figure;
%plot (srccat2_azi*180/pi, srccat2_el*180/pi);
%hold on
%plot (srcwsf2_azi*180/pi, srcwsf2_el*180/pi, 'r');
%set(gca, 'FontSize', 16);
%title ('Predicted Vs. observed track of CygA over 3 Hrs');
%xlabel ('Azimuth (deg)');
%ylabel ('Elevation (deg)');
%
%figure;
%plot (srccat3_azi(2:18)*180/pi, srccat3_el(2:18)*180/pi);
%hold on
%plot (srcwsf3_azi(2:18)*180/pi, srcwsf3_el(2:18)*180/pi, 'r');
%set(gca, 'FontSize', 16);
%title ('Predicted Vs. observed track of VirA over 3 Hrs');
%xlabel ('Azimuth (deg)');
%ylabel ('Elevation (deg)');

figure;
axis ([-3*180/pi 1.4*180/pi 0.3*180/pi 1*180/pi]);
plot (srccat1_azi*180/pi, srccat1_el*180/pi, 'LineWidth', 3);
hold on
plot (srcwsf1_azi*180/pi, srcwsf1_el*180/pi, 'r', 'LineWidth',3);
plot (srccat2_azi*180/pi, srccat2_el*180/pi);
plot (srcwsf2_azi*180/pi, srcwsf2_el*180/pi, 'r');
plot (srccat3_azi(2:18)*180/pi, srccat3_el(2:18)*180/pi);
plot (srcwsf3_azi(2:18)*180/pi, srcwsf3_el(2:18)*180/pi, 'r');
set(gca, 'FontSize', 16);
title ('Predicted Vs. observed track over 3Hrs, 21Sep2011, 1200UT');
xlabel ('Azimuth (deg)');
ylabel ('Elevation (deg)');

