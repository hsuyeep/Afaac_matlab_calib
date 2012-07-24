function C = khatrirao(A, B)

% C = khatrirao(A, B)
%
% Computes the Khatri-Rao (column wise Kronecker) product of A and B.
%
% Parameters
% A, B : 2-D matrices with equal second dimension. This is not checked by
%        the function, so entering matrices with insuitable dimensions may
%        produce unexpected results.
%
% Return value
% C    : column wise Kronecker product of A and B
%
% SJW, 2006

C = zeros(size(B, 1) * size(A, 1), size(A, 2));
for n = 1:size(A, 2)
    C(:, n) = kron(A(:, n), B(:, n));
end
