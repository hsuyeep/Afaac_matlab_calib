% Script to create primary beams for AARTFAAC images using SimElementBeam.m
% pep/12Oct15

obj = SimElementbeam (1, 'LBAOUTER', './', []);
lm = linspace (-1, 1, 1024);
freq = 45e6;
% for freq = [50e6:5e6:90e6]
   obj.createAvgbeam (1, lm, lm, [freq]); 
   obj.saveAARTFAACLMbeam (sprintf ('LBAOUTER_AARTFAAC_beamshape_%.02dMHz.hdf5',freq/1e6), 'hdf5');
% end;
