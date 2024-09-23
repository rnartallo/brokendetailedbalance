function [H,nH] = henrici(W)
%Calculate Henrici departure
n = norm(W,"fro");
[~,D] = eig(W);
eigenvalues = diag(D);
H = sqrt(n^2 - sum(abs(eigenvalues).^2));
nH = H/n;
end