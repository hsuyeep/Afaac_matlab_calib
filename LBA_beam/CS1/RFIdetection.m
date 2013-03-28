function clean = RFIdetection(acc, N)

nelem = size(acc, 1);
nch = size(acc, 3);
test = zeros(nch, 1);
treshold = sqrt(2 * nelem / N);

for ch = 1:nch
    Rhat = squeeze(acc(:, :, ch));
    Rwhite = Rhat ./ sqrt(diag(Rhat) * diag(Rhat).');
    test(ch) = norm(Rwhite, 'fro').^2;
end

clean = (abs([0; test; 0] - [test; 0; 0]) + abs([0; test; 0] - [0; 0; test])) < 10 * treshold;
