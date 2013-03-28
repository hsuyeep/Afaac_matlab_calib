% data for skymodel
load ../data/survey408MHz/intensity408MHz.mat
intensity408MHz = fliplr([intensity408MHz(:, ra >= 0), intensity408MHz(:, ra < 0)]);
intensity408MHz = conv2(intensity408MHz, ones(5), 'same');
ra = [ra(ra >= 0), ra(ra < 0) + 360 + 1e-8];
ra = ra * pi / 180;
dec = dec * pi / 180;

% other data and preparations
load cs10_fieldi.mat
freqidx = 11;
antidx = 1:2:96;
load CS10config_check20060920.mat
l = -1:0.01:1;
m = -1:0.01:1;
[lgrid, mgrid] = meshgrid(l, m);
dist = sqrt(lgrid.^2 + mgrid.^2);
thetai = asin(dist);
thetai(dist >= 1) = 0;
phii = mod(atan2(lgrid, mgrid), 2 * pi);
c = 2.99792e8;

% element beam pattern
for idxa = antidx
    idx1 = (idxa + 1) / 2;
    J_11a = FieldI(freqidx).E(:, :, 2, idxa);
    J_12a = FieldI(freqidx).E(:, :, 3, idxa);
    J_21a = FieldI(freqidx).E(:, :, 2, idxa+1);
    J_22a = FieldI(freqidx).E(:, :, 3, idxa+1);
    
    for idxb = antidx
        idx2 = (idxb + 1) / 2;

        J_11b = FieldI(freqidx).E(:, :, 2, idxb);
        J_12b = FieldI(freqidx).E(:, :, 3, idxb);
        J_21b = FieldI(freqidx).E(:, :, 2, idxb+1);
        J_22b = FieldI(freqidx).E(:, :, 3, idxb+1);
        
        E_11 = 0.5 * (J_11a .* conj(J_11b) + J_12a .* conj(J_12b));
        E_12 = 0.5 * (J_11a .* conj(J_21b) + J_12a .* conj(J_22b));
        E_21 = 0.5 * (J_21a .* conj(J_11b) + J_22a .* conj(J_12b));
        E_22 = 0.5 * (J_21a .* conj(J_21b) + J_22a .* conj(J_22b));
        
        lab = sin(FieldI(freqidx).THETA) .* cos(FieldI(freqidx).PHI);
        mab = sin(FieldI(freqidx).THETA) .* sin(FieldI(freqidx).PHI);
        delay = exp(-2 * pi * i * FieldI(freqidx).Freq / c * (lab * (xpos(CS10idx(idx1)) - xpos(CS10idx(idx2))) + mab * (ypos(CS10idx(idx1)) - ypos(CS10idx(idx2)))));
        if (idxb >= idxa)
            Eall_11(idx1, idx2, :, :) = E_11 .* delay;% * (idxb == idxa);
            Eall_12(idx1, idx2, :, :) = E_12 .* delay;% * (idxb == idxa);
            Eall_21(idx1, idx2, :, :) = E_21 .* delay;% * (idxb == idxa);
            Eall_22(idx1, idx2, :, :) = E_22 .* delay;% * (idxb == idxa);
        else
            Eall_11(idx1, idx2, :, :) = E_11 .* delay;
            Eall_12(idx1, idx2, :, :) = E_21 .* delay;
            Eall_21(idx1, idx2, :, :) = E_12 .* delay;
            Eall_22(idx1, idx2, :, :) = E_22 .* delay;
        end
    end
end
Eavg_11 = squeeze(mean(mean(Eall_11, 1), 2));
Eavg_12 = squeeze(mean(mean(Eall_12, 1), 2));
Eavg_21 = squeeze(mean(mean(Eall_21, 1), 2));
Eavg_22 = squeeze(mean(mean(Eall_22, 1), 2));
elembeam_xx = interp2(FieldI(freqidx).PHI, FieldI(freqidx).THETA, Eavg_11, phii, thetai);
elembeam_yy = interp2(FieldI(freqidx).PHI, FieldI(freqidx).THETA, Eavg_22, phii, thetai);

mask = 1.0 * (dist < 1);
weight = sqrt(1 - dist.^2);
weight(weight == 0) = Inf;

tstart = datenum(2007, 2, 20, 17, 06, 38);
dt = 0.01;
tstop = tstart + 53560 / 86400;
tobs = t_obs; %tstart:dt:tstop;

for idx = 1:length(tobs)
    disp([num2str((idx / length(tobs)) * 100) '% completed']);
    [alpha, delta] = lmtoradec(l, m, JulianDay(tobs(idx)));
    alpha(isnan(alpha)) = 0;
    delta(isnan(delta)) = 0;
    skyobs = interp2(ra, dec, intensity408MHz, alpha, delta);
    skyobs(isnan(skyobs)) = 0;
    [l0, m0] = radectolm(6.1237, 1.0267, JulianDay(tobs(idx))); % Cas A
    arraybeam = abs(xytolm(xpos(CS10idx(1:48)), ypos(CS10idx(1:48)), eye(48), l, m, c / FieldI(freqidx).Freq, l0, m0)).^2;
    beampowerx(idx) = sum(sum(skyobs .* elembeam_xx .* arraybeam .* mask .* (1 ./ weight)));
    beampowery(idx) = sum(sum(skyobs .* elembeam_yy .* arraybeam .* mask .* (1 ./ weight)));
end
