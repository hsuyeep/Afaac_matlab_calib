% Function to generate a complex cross-correlation matrix from values
% extracted from an MS using casa/listvis, and cleaned up in vim
% Extended to cater to all 4 polarizations, used with the superterp data.
% and to generate a .mat file
% pep/10Jan12
function [acm] = fill_acm_sterp (filename, nele)
fid = fopen (filename, 'r');
%fid2 = fopen ('fill_acm.deb', 'w');

% Read in the data
%% acc = textscan (fid, '%s %d %f %f %d %f %f %f');
ccm = textscan (fid, '%s %d %f %f %d %f %f %d %f %f %d %f %f %d %f %f %f');
ccm_amp = ccm{3};
ccm_ph = ccm{4};
%fprintf (fid2, 'Printing amps...\n');
%fprintf (fid2, '%f\n', acc_amp);
%fprintf (fid2, 'Printing phases...\n');
%fprintf (fid2, '%f\n', acc_ph);

%% load the vectored amplitude and phase into a matrix
begind = 1; endind = nele;
acc_amp_mat = zeros (nele);
acc_ph_mat = acc_amp_mat;

%%
for row = 1:nele
    %row
    %begind
    %endind
    acc_amp_mat (row, row:nele) = ccm_amp (begind:endind);
    acc_ph_mat  (row, row:nele) = ccm_ph (begind:endind);
    begind = endind + 1; endind = begind + nele - row - 1;
end

%%
acm = complex (acc_amp_mat.*cosd(acc_ph_mat), acc_amp_mat.*sind(acc_ph_mat));
newfilename = strrep (filename, '.txt', '.mat');
save (newfilename, 'acm');
%fprintf (fid2, 'Printing real part...\n');
%fprintf (fid2, '%f\n', real(acm(1, 1:48)));
fclose (fid);
%fclose (fid2);
     
