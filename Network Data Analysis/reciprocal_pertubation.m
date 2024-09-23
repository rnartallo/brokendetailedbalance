function [adj_eps] = reciprocal_pertubation(adj,eps)
N = length(adj);
adj_eps = adj;

for i =1:N
    for j = 1:N
        if adj(i,j)>0 && adj(j,i)==0
            adj_eps(j,i) = eps;
        end
    end
end
end