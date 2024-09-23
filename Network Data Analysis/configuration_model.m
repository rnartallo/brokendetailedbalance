function [B] = configuration_model(A)
out_degrees = sum((A>0), 2);
in_degrees = sum((A>0), 1);

[N,~] = size(A);
B = zeros(N);

inStubs = [];
outStubs = [];
for i = 1:N
    inStubs = [inStubs, repmat(i, 1, in_degrees(i))];
    outStubs = [outStubs, repmat(i, 1, out_degrees(i))];
end

inStubs = inStubs(randperm(length(inStubs)));
outStubs = outStubs(randperm(length(outStubs)));

B = zeros(N);

for i = 1:length(outStubs)
    u = outStubs(i);
    v = inStubs(i);
    B(u, v) = 1;
end