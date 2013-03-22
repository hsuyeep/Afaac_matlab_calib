function test_wracm2bin ()
	acc_wr = complex (randn (288), randn (288));
	% Make matrix hermitean symmetric, with reals on main diagonal.
	acc_wr(eye(288) == 1) = abs(diag(acc_wr)); 
	for ind=1:288
		for jind=ind+1:288
			acc_wr(ind, jind) = conj(acc_wr(jind, ind));
		end;
	end;
	tobs_wr = 12345678.23;
	freq_wr = 98765432.12;	
	flagants = [23, 56, 87, 230];

	fid = fopen ('/tmp/tmp_wracm2bintst.bin', 'wb');
	wracm2bin (fid, acc_wr, flagants, tobs_wr, freq_wr);
	fclose (fid);

	fid = fopen ('/tmp/tmp_wracm2bintst.bin', 'rb');
	
	% We always read in 288 elements anyway.
	[acc_rd, tobs_rd, freq_rd] = readms2float (fid, -1, -1, 288);
	% assertElementsAlmostEqual (acc_rd, acc_wr);
	% NOTE: assertElementAlmostEqual does not work as values are written out
	% and read in as floats.
	if (sum(sum(abs(real(acc_wr)) - abs(real(acc_rd)))) > 1e-3)
		error ('test_wracm2bin: acc real parts have significant error!');
	end;
	if (sum(sum(abs(imag(acc_wr)) - abs(imag(acc_rd)))) > 1e-3)
		error ('test_wracm2bin: acc imag parts have significant error!');
	end;
	assertEqual (tobs_rd, tobs_wr);
	assertEqual (freq_rd, freq_wr);
