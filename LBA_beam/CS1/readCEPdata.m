function data = readCEPdata(filename, nch, nstation, nsample, nsubband)

fid = fopen(filename, 'r');
rawdata = fread(fid, 'float');
fclose(fid);
npol = 4;

datare = reshape(rawdata(1:2:end), npol, nch, 0.5 * nstation * (nstation + 1), nsample, nsubband);
dataim = reshape(rawdata(2:2:end), npol, nch, 0.5 * nstation * (nstation + 1), nsample, nsubband);
size(datare)
for sample = 1:nsample
    for ch = 1:nch
        for sb = 1:nsubband
            for idx1 = 1:nstation
                for idx2 = idx1:nstation
                    baseline = (idx1-1)*nstation - 0.5*idx1*(idx1-1) + idx2
                    cohre = reshape(datare(:, ch, baseline, sample, sb), [2 2]).';
                    cohim = reshape(dataim(:, ch, baseline, sample, sb), [2 2]).';
                    acm(2*(idx2-1)+1:2*idx2, 2*(idx1-1)+1:2*idx1, (sb - 1) * nch + ch, sample) = (cohre + i * cohim);
                    if idx2 ~= idx1
                        acm(2*(idx1-1)+1:2*idx1, 2*(idx2-1)+1:2*idx2, (sb - 1) * nch + ch, sample) = (cohre + i * cohim)';
                    end
                end
            end
        end
    end
end
data = acm;