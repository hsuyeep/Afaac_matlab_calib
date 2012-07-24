% script to readin the VLSS catalog into a structure similar to the one
% used by SJWfor representing the 3CR catalog.
% pep/20Jan12


srcnumlimit = 100;% Only the first 'srcnumlimit' sources are picked
srcfluxlimit = 20; % In Jy. only sources having a flux above this are picked.
%['Could not find VLSS text catalog' fvlss] = fopen ('VLSScatalog07Jun26', 'r'); 
fvlss = fopen ('VLSScatalog07Jun26', 'r'); 

% Eliminate text in header, but print out relevant portions
str = fgets (fvlss)
while strcmp (strtok(str), 'RA(2000)') == 0
  str = fgets (fvlss);
  disp (str); % printf ('%s', str);
end

% eliminate  text line containing units
str = fgets (fvlss);

srccnt = 1;
while (feof (fvlss) == 0 && srccnt < srcnumlimit)
  str = fgets (fvlss);
  [val cnt] = sscanf (str, '%d %d %f %d %d %f %f', 7);
  
  if val(7) > srcfluxlimit
    srclistVLSS(srccnt).name = 'test';
    srclistVLSS(srccnt).alpha = val (1) + val (2)/60 + val (3)/3600;
    if val (4) < 0
      srclistVLSS(srccnt).delta = val (4) - val (5)/60 - val (6)/3600;
    else
      srclistVLSS(srccnt).delta = val (4) + val (5)/60 + val (6)/3600;
    end 
    srclistVLSS(srccnt).flux = val (7);
    srccnt = srccnt + 1;
  end 

  str = fgets (fvlss);    % ignore next line which only defines errors
end
save ('srclistVLSS.mat', 'srclistVLSS');
