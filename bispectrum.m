% Program to estimate the bispectrum of visibility timeseries.
% pep/06Feb12

clear all;

%%  Name of directory holding snapshot data matrices
dirname = 'sterp_tseries/';
% [~, lsout] = dos(['ls -1 ' dirname '*clean.txt']);
[~, lsout] = dos(['ls -1 ' dirname '*clean.mat']);

datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
Nelem = 288; % NOTE Hardcoded
Nbispect = 0.5 * (Nelem-1) * (Nelem-2);
bispect_amp = zeros (nfiles, Nbispect);
bispect_ph = bispect_amp;

%%
for fnum = 1:nfiles
    
    % read in the cross correlation matrix
    disp(['Processing file number ' num2str(fnum) ' of ' num2str(nfiles) '(' datafiles{1}{fnum} ')']);
    datestr (now)
    %acm = fill_acm_sterp (datafiles {1}{file}, 288);


    % Extract out the timestamp of the observation and load cross correlation
    % matrix
  
    ccm_fname = datafiles {1}{fnum};
    load (ccm_fname);

    [tok ccm_fname] = strtok (ccm_fname, '/');
    a = size (ccm_fname);
    while a(2) ~= 0
        [tok ccm_fname] = strtok (ccm_fname, '/');
        a = size (ccm_fname);
    end
    
    f_tobs = sscanf (tok, 'SB000_%2d%2d%2d')

    % the array covariance matrix (read in as the variable acc) is stored in upper triangular form, so we
    % complete it
    acc = acm + acm' - diag(diag(acm));
    Nelem = size(acc, 1);
    
    %% compute the bispectrum only for unique antenna triples, for all unique sets of 3 antennas
    cnt = 0;
    %Nelem = 6;
    
    flag = zeros (Nelem);
    for i = 1:Nelem
        for j = i+1:Nelem
            for k=j+1:Nelem
                %  disp ([num2str(i) ',' num2str(j) '/' num2str(j) ','
                %  num2str(k) '/' num2str(k) ',' num2str(i)]);
                if (flag (i, j) == 0 || flag (j, k) == 0 || flag (k, i) == 0)
                    flag (i,j) = 1; flag (j,i) = 1; 
                    flag (j,k) = 1; flag (k,j) = 1;
                    flag (k,i) = 1; flag (i,k) = 1;
                    cnt = cnt + 1;
                    bispect_amp (fnum, cnt) = abs  (acc(i,j)) * abs  (acc(j,k)) * abs  (acc(k,i));
                    bispect_ph  (fnum, cnt) = angle(acc(i,j)) + angle(acc(j,k)) + angle(acc(k,i));
                   
                end
                
            end
        end
    end
    cnt
    

end

