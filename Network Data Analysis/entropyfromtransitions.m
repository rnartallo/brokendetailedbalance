function [Phi] = entropyfromtransitions(P)
[N,~] = size(P);
M = 0.5*(P+P');
Phi =0;
for i=1:N
    for j=1:N
        if (P(i,j)>0)
            Phi = Phi + 0.5*P(i,j)*log(P(i,j)/M(i,j));
        end
        if (P(j,i)>0)
            Phi = Phi + 0.5*P(j,i)*log(P(j,i)/M(i,j));
        end
    end
end
end