clear all

dirname = '../../data/DE604C/xcstat64s_change3/';
nelem = 32;
nyquist = 100;

[dummy, lsout] = dos(['ls -1 ' dirname]);
datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
offset = length(dirname);

acc = zeros(nelem, nelem, nfiles);
for idx = 1:nfiles
    disp(['reading file ' datafiles{1}{idx} '(' num2str(idx) ' of ' num2str(nfiles) ')']);
    fid = fopen([dirname datafiles{1}{idx}], 'r');
    data = fread(fid, 'double');
    fclose(fid);
    cdata = data(1:2:2*nelem^2) + i * data(2:2:2*nelem^2);
    accadd = reshape(cdata, [nelem, nelem]);
    acc(:, :, idx) = accadd;% + 0 * accadd' .* (accadd == 0);
end
%freq = 0:nyquist/nfiles:nyquist;
%freq = freq(1:end-1) * 1e6;
