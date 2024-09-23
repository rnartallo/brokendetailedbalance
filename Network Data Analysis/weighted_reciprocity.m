function [r] = weighted_reciprocity(W)
[N,~] = size(W);
W_symm = min(W,W');
num = sum(sum(W_symm-diag(diag(W_symm))));
denom = sum(sum(W))- sum(diag(W));
r = num/denom;
end