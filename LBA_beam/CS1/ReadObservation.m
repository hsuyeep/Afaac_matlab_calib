%!/usr/bin/octave

InputFile = 'vis.dat';
NoSubbands = 2;
BaseLines = 3; % including autocorrelations
TimeSteps = 1;
Polarisations = 4;
ChannelsPerSubband = 1; % is always 1 for now

% the 2 is because of the 2 doubles per complex<double>
datashape = [2, Polarisations, ChannelsPerSubband, BaseLines, TimeSteps, NoSubbands]
NoVis = TimeSteps * NoSubbands * BaseLines * ChannelsPerSubband * Polarisations

fid = fopen(InputFile, 'r');
rawdata = fread(fid, NoVis * 2, 'float');
size(rawdata)


data = reshape (rawdata, datashape);

c = 1 % there is only 1 channel
for t = 1:TimeSteps
  for s = 1:NoSubbands
    for b = 1:BaseLines
      disp(sprintf('TimeStep %6d - Subband %2d - BaseLine %1d    ', t, s, b));
      disp(sprintf('POLS (XX, XY, YX, YY): '));
      for p = 1:Polarisations
        disp(sprintf('(%f, %f)', data(1, p, c, b, t, s), data(2, p, c, b, t, s)));
      end
      disp(sprintf('\n'));
    end
  end
end

