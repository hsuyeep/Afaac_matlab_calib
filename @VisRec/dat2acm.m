% Function to convert data stored in the xx/yy linear arrays into
% a square symmetric Array Covariance Matrix.
function obj = dat2acm (obj)
    assert (~isempty (obj.xx));
    assert (~isempty (obj.yy));
    t = triu (ones(obj.nelem));
    obj.acm_xx = zeros (obj.nelem);
    obj.acm_xx (t == 1) = obj.xx; 
    acm_diag = diag (diag (obj.acm_xx)); 
    obj.acm_xx = conj (obj.acm_xx + obj.acm_xx' - acm_diag);

    obj.acm_yy = zeros (obj.nelem);
    obj.acm_yy (t == 1) = obj.yy; 
    acm_diag = diag (diag (obj.acm_yy)); 
    obj.acm_yy = conj (obj.acm_yy + obj.acm_yy' - acm_diag);
end;
