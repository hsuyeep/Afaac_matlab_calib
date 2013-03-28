% requires cs10_fieldi.mat

% parameters
component = 3;
antenna = 1:2:96;
freqidx = 11;

l = sin(FieldI(freqidx).THETA) .* cos(FieldI(freqidx).PHI);
m = sin(FieldI(freqidx).THETA) .* sin(FieldI(freqidx).PHI);
c = 2.99798e8;

disp(['Frequency: ' num2str(FieldI(freqidx).Freq / 1e6) ' MHz']);

for idx = antenna
    idx1 = (idx + 1) / 2;
    disp(['position:  ' num2str(xpos(CS10idx(idx1))) ' ' num2str(ypos(CS10idx(idx1)))]);
    disp(['distance:  ' num2str(sqrt((xpos(CS10idx(idx1)) - xpos(CS10idx(1)))^2 + (ypos(CS10idx(idx1)) - ypos(CS10idx(1)))^2))]);
    disp(['direction: ' num2str(atan2(xpos(CS10idx(idx1)) - xpos(CS10idx(1)), ypos(CS10idx(idx1)) - ypos(CS10idx(1))) * 180 / pi)]);
    
    delay = exp(-2 * pi * i * (FieldI(freqidx).Freq / c) * ((xpos(CS10idx(idx1)) - xpos(CS10idx(1))) * l + (ypos(CS10idx(idx1)) - ypos(CS10idx(1))) * m));
    
    J_11 = FieldI(freqidx).E(:, :, 2, idx);
    J_12 = FieldI(freqidx).E(:, :, 3, idx);
    J_21 = FieldI(freqidx).E(:, :, 2, idx+1);
    J_22 = FieldI(freqidx).E(:, :, 3, idx+1);
    
    J1_11 = FieldI(freqidx).E(:, :, 2, 1);
    J1_12 = FieldI(freqidx).E(:, :, 3, 1);
    J1_21 = FieldI(freqidx).E(:, :, 2, 2);
    J1_22 = FieldI(freqidx).E(:, :, 3, 2);
    
    E_11 = 0.5 * (J_11 .* conj(J1_11) + J_12 .* conj(J1_12)) .* delay;
    E_12 = 0.5 * (J_11 .* conj(J1_21) + J_12 .* conj(J1_22)) .* delay;
    E_21 = 0.5 * (J_21 .* conj(J1_11) + J_22 .* conj(J1_12)) .* delay;
    E_22 = 0.5 * (J_21 .* conj(J1_21) + J_22 .* conj(J1_22)) .* delay;
    
    E1_11 = 0.5 * (J1_11 .* conj(J1_11) + J1_12 .* conj(J1_12));
    E1_12 = 0.5 * (J1_11 .* conj(J1_21) + J1_12 .* conj(J1_22));
    E1_21 = 0.5 * (J1_21 .* conj(J1_11) + J1_22 .* conj(J1_12));
    E1_22 = 0.5 * (J1_21 .* conj(J1_21) + J1_22 .* conj(J1_22));
    
    subplot(2, 2, 1);
    pcolor(l, m, abs(E_11 - E1_11) / max(abs(E1_11(:))) );
    shading flat
    colorbar
    
    subplot(2, 2, 2);
    pcolor(l, m, abs(E_12 - E1_12) / max(abs(E1_12(:))) );
    shading flat
    colorbar
    
    subplot(2, 2, 3);
    pcolor(l, m, abs(E_21 - E1_21) / max(abs(E1_21(:))) );
    shading flat
    colorbar
    
    subplot(2, 2, 4);
    pcolor(l, m, abs(E_22 - E1_22) / max(abs(E1_22(:))) );
    shading flat
    colorbar
    
    Eall_11(idx1, :, :) = E_11;
    Eall_12(idx1, :, :) = E_12;
    Eall_21(idx1, :, :) = E_21;
    Eall_22(idx1, :, :) = E_22;
    
    pause(0.1)
end

