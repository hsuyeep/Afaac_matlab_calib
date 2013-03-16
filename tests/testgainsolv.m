% Script to test the gainsolv.m function, either with an observed dataset, or
% with a more detailed simulated sky.
% Pep/21Dec12

function testgainsolv ()

% Simulated instantaneous complex voltage from each antenna.
ants = [10+20i 14+25i 12+22i]; 

nants = length(ants);
model = ants' * ants; % Model visibilities
model = model - diag(diag(model));

% Artificially introduced gains.
gains = [10+15i 24+10i 42+12i]
data  = diag(gains)*model*diag(gains)';
recov = gainsolv (1e-6, model, data, ones (nants, 1))

disp ('Original gain amp:')
abs(gains' * gains) 
disp ('Original gain phase:')
angle (gains' * gains)

disp ('Recovered gain amp:')
abs(recov  * recov') 
disp ('Recovered gain phase:')
angle (recov  * recov')

