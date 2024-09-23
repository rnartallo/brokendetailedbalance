function [F,h] = trophic_coherence(W)
[N,~] = size(W);
col_sum = sum(W); row_sum = sum(W');
u = col_sum+row_sum;
v = col_sum-row_sum;
Lambda = diag(u)-W-W';

% modify system to have a unique solution

% % calculate traditional laplacian 
% Laplace = diag(sum((W>0),2)) - (W>0);
% % nullity is number of connected components; basis vectors are constant on
% % CC
% null_basis =null(Laplace);
% [~,no_cc] = size(null_basis);
% % use kmeans (without actually clustering - using exact number of states) to sort the vector into
% % connected components
% %connected_components = kmeans(null_basis(:,1),no_cc);
connected_components = conncomp(graph(0.5*(W+W')));
no_cc = length(unique(connected_components));

for n=1:no_cc
        % find first node in a connected component and modify its equation
        firsteq_cc = find(connected_components==n,1);
        Lambda(firsteq_cc,:) = zeros(N,1);
        Lambda(firsteq_cc,firsteq_cc) = 1;
        v(firsteq_cc) = 1;
end

% solve system for unique solution
h = Lambda\v';

% normalise lowest level in each cc to be 0
for n=1:no_cc
    h(connected_components==n) = h(connected_components==n)-min(h(connected_components==n));
end

num = 0;
denom =sum(W-diag(diag(W)),'all');
for i=1:N
    for j=1:N
        if i~=j
            num = num + W(i,j)*((h(j)-h(i)-1)^2);
        end
    end
end
F = num/denom;
end