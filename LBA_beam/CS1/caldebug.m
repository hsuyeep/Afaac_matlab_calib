% requires session20050823-1.mat
load srclist3C

for fidx = 1:length(freq)
    if (clean(fidx))
        disp(['calibrating subband ' num2str(fidx)]);
    [flux, sortidx] = sort(srclist3C.flux, 'descend');
    alpha = srclist3C.alpha(sortidx);
    delta = srclist3C.delta(sortidx);
    names = srclist3C.names(sortidx);
    alpha = alpha(1:4);
    delta = delta(1:4);
    names = names(1:4);
    flux = flux(1:4);

    [l, m] = radectolm(alpha, delta, JulianDay(datenum(2005, 3, 23, 13, 25, 23.16)));

    k = 2 * pi * freq(fidx) / 2.99792e8;

    load itsbeam.mat
    lgrid = -1:0.04:1;
    mgrid = -1:0.04:1;
    [li, mi] = meshgrid(lgrid, mgrid);
    phi = atan2(li, mi) * 180 / pi;
    phi(phi < 0) = phi(phi < 0) + 360;
    theta = asin(sqrt(li.^2 + mi.^2)) * 180 / pi;
    for idx = 1:length(Freq);
        elembeam(:, :, idx) = interp2(Phi, Theta, DirPhi(:, :, idx) + DirTheta(:, :, idx), phi, real(theta));
    end

    R0 = zeros(length(xpos), length(xpos));
    for idx = 1:4
        if ~isnan(l(idx))
            asrc = exp(-i * k * (l(idx) * xpos + m(idx) * ypos)).';
            freqval = [20 30 50 80 120 180 240] * 1e6;
            fluxval = (178e6 ./ freqval).^0.7 * flux(idx);
            sigma_src = interp1(freqval, fluxval, freq(fidx), 'linear');
            att = interp3(lgrid, mgrid, Freq, elembeam, l(idx), m(idx), freq(fidx)/1e9, 'linear');
            R0 = R0 + att * sigma_src * asrc * asrc';
        end
    end
    refacm(fidx, :, :) = R0;
    u = meshgrid(xpos) - meshgrid(xpos).';
    v = meshgrid(ypos) - meshgrid(ypos).';
    mask = sqrt(u.^2 + v.^2) > 40;

    alpha = computeAlpha(acm(:, :, fidx) .* mask, R0 .* mask);
    [v, d] = eig(alpha);
    [dummy, idx] = max(max(d));
    gain = v(:, idx);

    Rtest = diag(gain) * (R0 .* mask) * diag(gain)';
    compidx = (Rtest ~= 0);
    acmfidx = acm(:, :, fidx);
    norm = 1 / sqrt(mean(abs(Rtest(compidx) ./ acmfidx(compidx))));
    gain = norm * gain / (gain(1) / abs(gain(1)));
    result(fidx, :) = gain;
    else
        disp(['skipping subband ', num2str(fidx)]);
        result(fidx, :) = zeros(1, length(xpos));
    end
end
