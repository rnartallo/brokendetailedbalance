% Load the network you want to analyse

[henrici_original,normalised_henrici_original] = henrici(adj);
irreciprocity_original = 1-weighted_reciprocity(adj);
[F,~] = trophic_coherence(adj);
trophic_directedness_original=sqrt(1-F);

Theta = 0.5;
gamma = 0.0025;
rwdt = 100;
sigma = 0.75;
iterations = 100;

for i =1:100
    adj_null = configuration_model(adj);
    [h,nh]= henrici(adj_null);
    henrici_null(i) = h;
    normalised_henrici_null(i)=nh;
    irreciprocity_null(i) = 1- weighted_reciprocity(adj_null);
    [F,~] = trophic_coherence(adj_null);
    trophic_directedness_null(i)=sqrt(1-F);
    
    % OU EPR
    I = eye(N);
    B = Theta*(I-gamma*adj_null);
    [~,D]=eig(B);
    flag = all(real(diag(D))>0);
    if flag
        S = lyap(B,-2*sigma*I);
        Q = B*S-sigma*I;
        OUP(i) = -trace((1/sigma)*B*Q);
    else
        i
    end
    
    % Ising EPR
    Phi = zeros(iterations);
    m_prev = ones([N,1]);
    for t =1:iterations
        v = adj*m_prev;
        m_new = tanh(v);
        for k=1:N
            for j=1:N
                D(k,j)=adj_null(k,j)*(1-m_new(k)^2)*(1-m_prev(j)^2);
                Phi(t) = Phi(t)+ (adj_null(k,j)-adj_null(j,k))*D(k,j);
            end
        end
        m_prev = m_new;
    end
    Ising(i)=Phi(iterations);

    cs = conncomp(digraph(adj_null),'Type','strong');
    cw = conncomp(digraph(adj_null),'Type','weak');

    if length(unique(cs))~=1 && length(unique(cw))==1
        adj_null = reciprocal_pertubation(adj_null,0.01);
    end

    if length(unique(cw))~=1
        i
    end
    
    s = sum(adj_null,2);
    for k =1:length(adj_null)
        for j=1:length(adj_null)
            T(k,j) = adj_null(k,j)/s(k);
        end
    end
    [vec,val]=eig(T');
    ss = abs(vec(:,1));
    P = rwjointtransition(ss,T);
    DTRW(i) = entropyfromtransitions(P);
end

%%
writematrix(henrici_null,'LondonHenriciNull.txt')
writematrix(normalised_henrici_null,'LondonNormalisedHenriciNull.txt')
writematrix(irreciprocity_null,'LondonIrreciprocityNull.txt')
writematrix(trophic_directedness_null,'LondonTrophicDNull.txt')

writematrix(OUP,'LondonOUNull.txt')
writematrix(Ising,'LondonIsingNull.txt')
writematrix(DTRW,'LondonRWNull.txt')


