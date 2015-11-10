% Script to generate variation of noise as a function of RA/DEC
% using the 24hr observation data.
% pep/04Nov15

function [tobs, noisets] = func_24hrnoisemap ()

    % List of filenames containing the image data.
    fnames = {'~/WORK/AARTFAAC/Reobs/20Nov13/r01/SB002_LBA_OUTER_8b2sbr01_1ch_1_convcal_1_noel_fftimg.bin.tmpmove', '~/WORK/AARTFAAC/Reobs/20Nov13/r02/SB002_LBA_OUTER_8b2sbr02_1ch_7560_convcal_1_noel_fftimg.bin', '~/WORK/AARTFAAC/Reobs/20Nov13/r02/SB002_LBA_OUTER_8b2sbr02_1ch_9860_convcal_1_noel_fftimg.bin', '~/WORK/AARTFAAC/Reobs/20Nov13/r03/SB002_LBA_OUTER_8b2sbr03_1ch_1337-1850_3859_convcal_1_noel_fftimg.bin', '~/WORK/AARTFAAC/Reobs/20Nov13/r04/SB002_LBA_OUTER_8b2sbr04_1ch_0006-0442_1_convcal_1_noel_fftimg.bin', '~/WORK/AARTFAAC/Reobs/20Nov13/r04/r04/SB002_LBA_OUTER_8b2sbr04_1ch_0006-0442_1_convcal_1_noel_fftimg.bin', '~/WORK/AARTFAAC/Reobs/20Nov13/r04/r04/SB002_LBA_OUTER_8b2sbr04_1ch_0006-0442_8812_convcal_1_noel_fftimg.bin'}; % , '~/WORK/AARTFAAC/Reobs/20Nov13/r04/r04/};

    % Give out file information
    for ind = 1:length(fnames)
        fimg = fopen (fnames{ind}, 'r');
        img = readimg2bin (fimg, 0);
		recsize = 8+ ... % tobs
				  8+ ... % freq
				  4+ ... % pix2laxis
				  4*img.pix2laxis + ... % laxis
				  4+...  % pix2maxis
				  4*img.pix2maxis + ... % maxis
			      4*img.pix2laxis*img.pix2maxis; % img contents.
        
        fseek (fid, 'eof', -recsize);
        imglast = readimg2bin(fimg, 0);
        fprintf (2, '<-- File %s\n', fnames{ind});
        fprintf (2, '<-- %dx%d image, freq: %f\n', img.pix2laxis, img.pix2maxis, img.fobs);
        fprintf (2, '<-- Timerange : %s - %s', datestr (mjdsec2datenum(img.tobs)), datestr (mjdsec2datenum(imglast.tobs));
    end;

    % Initially, we look at the noise behaviour only at the zenith
    for ind = 1:length (fnames)
        
    end;
    

    

