function [M,N] = Conv_struc2matix(S,m,p)

% Converts a struct of vercotrs of indexes into a matrix of logicals

if nargin == 1
    % find the max index
    a = [];
    for i = 1:length(S)
        a = [a S(i).val];
    end
    my = max(a);
    mx = length(S);
elseif nargin == 3
    mx = m;
    my = p;
end

M = sparse(mx,my);
for i = 1:length(S)
    M(i,S(i).val) = 1;
    N(i) = length(S(i).val);
end
% nnz(M)

