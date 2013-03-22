%
function test_readwrcalsol2bin ()
	% Open test data.
	testfile = 'unittest_SB004_LBA_OUTER.bin';
	fid = fopen (testfile, 'rb');
	fsol = fopen ('/tmp/testsol.bin', 'wb');

	try
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	catch ME
		error ('unittest: wrcalsol2bin.m FAILED!\n');
	end;

	% Generate cal solution
	currflagant = [51, 238, 273]; 
	try
		sol = pelican_sunAteamsub (acc, tobs, freq, eye(length(acc)), ... 
							currflagant, 0, 1, [], []);
	catch ME
		error ('unittest: pelican_sunAteamsub.m FAILED!\n');
	end;
	
	% write out solution.
	try
		wrcalsol2bin (fsol, sol);
	catch ME
		error ('unittest: wrcalsol2bin.m FAILED!\n');
	end;
	fclose (fsol);
	

	fsol = fopen ('/tmp/testsol.bin', 'rb');
	% Read it back in.
	try
		rdsol = readcalsol (fsol);
	catch ME
		error ('unittest: readcalsol.m FAILED!\n');
	end;
	fclose (fsol);
	
	% assertions
	assertEqual (sol.tobs, rdsol.tobs);
	assertEqual (sol.freq, rdsol.freq);
	assertEqual (single (sol.flagant), rdsol.flagant');
	assertElementsAlmostEqual (single(real(sol.gainsol)), rdsol.real_gainsol');
	assertElementsAlmostEqual (single(imag(sol.gainsol)), rdsol.imag_gainsol');
	assertElementsAlmostEqual (sol.sigmas, rdsol.sigmas);
	assertElementsAlmostEqual (sol.thsrc_wsf, rdsol.thsrc_wsf);
	assertElementsAlmostEqual (sol.phisrc_wsf, rdsol.phisrc_wsf);
	assertElementsAlmostEqual (single (real (sol.sigman(:))),rdsol.real_sigman);
	assertElementsAlmostEqual (single (imag (sol.sigman(:))),rdsol.imag_sigman);
	
