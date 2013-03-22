% Function to generate a complex cross-correlation matrix from values
% extracted from an MS using casa/listvis, and cleaned up in vim
% pep, Jan12

function [acm] = fill_acm (filename, nele)
fid = fopen (filename, 'r');
fid2 = fopen ('fill_acm.deb', 'w');

% Read in the data
acc = textscan (fid, '%s %d %f %f %d %f %f %f');
acc_amp = acc{3};
acc_ph = acc{4};
fprintf (fid2, 'Printing amps...\n');
fprintf (fid2, '%f\n', acc_amp);
fprintf (fid2, 'Printing phases...\n');
fprintf (fid2, '%f\n', acc_ph);

% load the vectored amplitude and phase into a matrix
begind = 1; endind = nele;
acc_amp_mat = zeros (nele+1);
acc_ph_mat = acc_amp_mat;

for i = 1:nele
     acc_amp_mat (i, i+1:nele+1) = acc_amp (begind:endind);
     acc_ph_mat  (i, i+1:nele+1) = acc_ph (begind:endind);
     begind = endind + 1; endind = begind + nele - i - 1;
end

acm = complex (acc_amp_mat.*cosd(acc_ph_mat), acc_amp_mat.*sind(acc_ph_mat));
fprintf (fid2, 'Printing real part...\n');
fprintf (fid2, '%f\n', real(acm(1, 1:48)));
fclose (fid);
fclose (fid2);
     
