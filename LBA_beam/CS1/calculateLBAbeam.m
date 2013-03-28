function elembeamx = calculateLBAbeam(l, m, freq, antno)

load cs10_fieldi.mat
[lgrid, mgrid] = meshgrid(l, m);
dist = sqrt(lgrid.^2 + mgrid.^2);
thetai = asin(dist);
thetai(dist >= 1) = 0;
phii = mod(atan2(lgrid, mgrid), 2 * pi);
E_11 = zeros(length(antno), size(FieldI(1).PHI, 1), size(FieldI(1).PHI, 2), length(FieldI(1).Freq));
simfreq = zeros(1, length(FieldI));
for fidx = 1:length(FieldI)
    for idx = 1:length(antno)
        J_11 = FieldI(fidx).E(:, :, 2, antno(idx));       
        J_12 = FieldI(fidx).E(:, :, 3, antno(idx));
        E_11(idx, :, :, fidx) = 0.5 * (J_11 .* conj(J_11) + J_12 .* conj(J_12));
    end
    simfreq(fidx) = FieldI(fidx).Freq;
end
Eavg11 = squeeze(mean(E_11, 1));
elembeamx = zeros(length(l), length(m), length(freq));
for idx = 1:length(freq)
    elembeamx(:, :, idx) = interp3(FieldI(1).PHI(1, :), FieldI(1).THETA(:, 1), simfreq, Eavg11, phii, thetai, freq(idx) * ones(size(phii))) .* (dist < 1);
end
