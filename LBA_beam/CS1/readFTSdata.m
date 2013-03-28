clear all

dirname = './';
nelem = 96;
nch = 512;
nyquist = 80;

[dummy, lsout] = dos(['ls -1 ' dirname '*acc*.dat']);
datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
offset = length(dirname);

acc = zeros(nfiles, nelem, nelem, nch);
for idx = 1:nfiles
    disp(['reading file ' datafiles{1}{idx}]);
    fid = fopen(datafiles{1}{idx}, 'r');
    data = fread(fid, 'double');
    fclose(fid);
    redata = data(1:2:end);
    imdata = data(2:2:end);
    snapshotacc = reshape(redata, [nelem, nelem, nch]) + i * reshape(imdata, [nelem, nelem, nch]);
    acc(idx, :, :, :) = snapshotacc(:, :, :);
    timestamp(idx) = datenum(datafiles{1}{idx}(offset+1:offset+15), 'yyyymmdd_HHMMSS');
end
freq = 0:nyquist/nch:nyquist;
freq = freq(1:end-1) * 1e6;
