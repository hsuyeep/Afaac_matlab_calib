% Script to test readms2float.m
% pep/21Mar13

% Location of test file. Should contain at least 10-20 timeslices due to the
% windowing of timeslices in certain functions.

function test_readms2float ()
	testfile = 'unittest_SB004_LBA_OUTER.bin';

	% - Basic reading and writing of .bin files.
	fid = fopen (testfile, 'rb');
	try
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	catch ME
		error ('unittest: readms2float.m FAILED!\n');
	end;
	
