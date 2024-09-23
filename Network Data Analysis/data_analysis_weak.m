myDir = uigetdir; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all mat files in struct

L = length(myFiles);
normalised_henrici = zeros(L,1);
henrici_departure = zeros(L,1);
directionality = zeros(L,1);
reciprocity = zeros(L,1);
network_size = zeros(L,1);
strongc = 0;
weakc = 0;

normalised_henrici = zeros(L,1);
henrici_departure = zeros(L,1);
directionality = zeros(L,1);
reciprocity = zeros(L,1);

for k = 1:length(myFiles)
  clear T
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  %fprintf(1, 'Now reading %s\n', fullFileName);
  load(baseFileName)
  cs = conncomp(digraph(adj),'Type','strong');
  cw = conncomp(digraph(adj),'Type','weak');
  N = length(adj);
  if length(unique(cs))==1
      fprintf(1, 'Strongly connected: ', fullFileName);
      fprintf(1, 'Now reading %s\n', fullFileName);
      k
      strong =true;
      strongc = strongc+1;
      [H,nH] = henrici(adj);
      normalised_henrici(k) = nH;
      henrici_departure(k) = H;
      [F,~] = trophic_coherence(adj);
      directionality(k) = 1-F;
      reciprocity(k) = weighted_reciprocity(adj);

      s = sum(adj,2);
      for i =1:length(adj)
          for j=1:length(adj)
              T(i,j) = adj(i,j)/s(i);
          end
      end

      [vec,val]=eig(T');
      ss = abs(vec(:,1));
      P = rwjointtransition(ss,T);
      DTRW(k) = entropyfromtransitions(P);

      sigma = 0;
      B = bidirectional_transition_matrix(N,T);
      for i =1:N
          for j = 1:N
              if B(i,j)>0
                sigma = sigma +  B(j,i)*ss(j)*log((B(j,i)*ss(j))/(B(i,j)*ss(i)));
              end
          end
      end
      CTRW(k) = sigma;

      

  elseif length(unique(cw))==1
      % Connected
      adj_eps = reciprocal_pertubation(adj,0.01);
      [H,nH] = henrici(adj_eps);
      normalised_henrici(k) = nH;
      henrici_departure(k) = H;
      [F,~] = trophic_coherence(adj_eps);
      directionality(k) = 1-F;
      reciprocity(k) = weighted_reciprocity(adj_eps);

      s = sum(adj_eps,2);
      for i =1:length(adj_eps)
          for j=1:length(adj_eps)
              T(i,j) = adj_eps(i,j)/s(i);
          end
      end

      [vec,val]=eig(T');
      ss = abs(vec(:,1));

      P = rwjointtransition(ss,T);
      DTRW(k) = entropyfromtransitions(P);

      sigma = 0;
      B = bidirectional_transition_matrix(N,T);
      for i =1:N
          for j = 1:N
              if B(i,j)~=0 && B(j,i)~=0
                sigma = sigma +  B(j,i)*ss(j)*log((B(j,i)*ss(j))/(B(i,j)*ss(i)));
              end
          end
      end
      CTRW(k) = sigma;
  else
      fullFileName
      CTRW(k) = nan;
      DTRW(k) = nan;
      [H,nH] = henrici(adj);
      normalised_henrici(k) = nH;
      henrici_departure(k) = H;
      [F,~] = trophic_coherence(adj);
      directionality(k) = 1-F;
      reciprocity(k) = weighted_reciprocity(adj);
  end

end


%%

