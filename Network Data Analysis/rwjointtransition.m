function [P] = rwjointtransition(pi,T)
N = length(pi);
P = zeros(N,N);
for i=1:N
    for j=1:N
        P(i,j) = T(i,j)*pi(i);
    end
end
end