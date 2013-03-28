clear skymap
chsel = 1:216;
skymap = zeros(length(chsel), 30, 30);
idx = 1;
offset = 29;
for ch = chsel
    measidx = offset;
    for mval = m;
        for lval = l;
            if ((lval^2 + mval^2) <= 1)
                skymap(ch-chsel(1)+1, find(m == mval), find(l == lval)) = bst(ch, measidx);
                measidx = measidx + 1;
            else
                skymap(ch-chsel(1)+1, find(m == mval), find(l == lval)) = NaN;
            end 
        end
    end
end
