% Please kindly cite the paper Junyi Guan, Sheng Li, Xiongxiong He, and Jiajia Chen,
%"Peak-graph-based fast density peak clustering for image segmentation,"
% IEEE SIGNAL PROCESSING LETTERS, 2021,Doi:?
% The code was written by Junyi Guan in 2020.
function [CL,Centers] = PGDPC(data,k)
close all
[n,d]  = size(data); % N: the number of data points; d:dimensions
%% Normalization
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;
%% Fast KNN calculation based on Kd-tree (when dimension is not large than 10)
if d<=11;[knn_pt,knn_dist] = knnsearch(data,data,'k',k);else;dist = squareform(pdist(data));[knn_dist,knn_pt] = sort(dist,2);end
%% FKNN-DPC Density Evaluation method based on KNN;
%% "Xie,Juanying, et al, "Robust clustering by detecting density peaks and assigning points based on fuzzy weighted K-nearest neighbors"
rho=sum(exp(-knn_dist(:,2:k)),2)';
%% Allocation
[~,ordrho]=sort(rho,'descend');
RNs = [];%RNs:root nodes
for i=1:n
    depth(ordrho(i))= 0;
    for j=2:k
        neigh=knn_pt(ordrho(i),j);
        if(rho(ordrho(i))<rho(neigh))
            PN(ordrho(i))=neigh; %PN:parent node
            depth(ordrho(i))= depth(neigh) + 1;
            break
        end
    end
    if depth(ordrho(i))==0
        RNs = [ordrho(i) RNs];
    end
end
%% label initialization
for i=1:n; subL(i)=-1; end %% subL:sub label from subtrees (subclusters)
%% sub label Assignation
NRNs = length(RNs);%% NRNs:Number of root nodes
subL(RNs) = (1:NRNs);
for i=1:n
    if (subL(ordrho(i))==-1)
        subL(ordrho(i))=subL(PN(ordrho(i)));
    end
end
%% Obiain edges matrix
rho_RNs= rho(RNs); [~,ordrho_RNs]=sort(rho_RNs,'descend');
edges = Inf*ones(NRNs,NRNs);% edges: Edge matrix between subtrees
for i=1:n
    for j = 2:k
        jj = knn_pt(i,j);
        if subL(i)~=subL(jj) & edges(subL(i),subL(jj))> depth(i)+depth(jj)+1
            if find(knn_pt(jj,2:k)==i);
                edges(subL(i),subL(jj)) = depth(i)+depth(jj)+1; edges(subL(jj),subL(i)) = depth(i)+depth(jj)+1;
            end
        end
    end
end
%% Construct Peak-Graph
G=sparse(NRNs,NRNs);
for i=1:NRNs
    for j = 1:NRNs
        if edges(i,j) ~= Inf G(i,j) = edges(i,j); G(j,i) = edges(i,j); end
    end
end
[d_dijks ~] = dijkstra(G,(1:NRNs));
d_dijks = round(d_dijks);
Delta_RNs = Inf*ones(NRNs,1);
PN_RNs  = -1*ones(NRNs,1);
for i = 2:NRNs
    ii = ordrho_RNs(i);
    for j = 1:i-1
        jj = ordrho_RNs(j);
        if Delta_RNs(ii) > d_dijks(ii,jj)
            Delta_RNs(ii) = d_dijks(ii,jj);
            PN_RNs(ii) = jj;
        end
    end
end
delta = zeros(n,1);
delta(RNs) =  Delta_RNs;
if NRNs > 1
    delta(delta==Inf) = max(delta(delta~=Inf))*1.5;
end
%% center selection
figure(1);
plot(rho,delta,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
title ('Decision Graph','FontSize',15.0);
xlabel ('\rho');
ylabel ('\delta');
rect = getrect;
Rhomin=rect(1);
Deltamin=rect(2);
NCLUST=0;
for i=1:NRNs
    CL_RNs(i)=-1;
end
for i=1:NRNs
    if rho_RNs(i)>Rhomin & Delta_RNs(i)>Deltamin
        NCLUST=NCLUST+1;
        CL_RNs(i)=NCLUST;
        icl(NCLUST)=i;
    end
end
%% Assignation
for i=1:NRNs
    if (CL_RNs(ordrho_RNs(i))==-1)
        CL_RNs(ordrho_RNs(i))=CL_RNs(PN_RNs(ordrho_RNs(i)));
    end
end
for i=1:NRNs
    CL(subL== i) = CL_RNs(i);
end
CL = CL';
Centers = RNs(icl);
%% show result
cmap = colormap;
for i=1:NCLUST
    ic=int8((i*64.)/(NCLUST*1.));
    plot(rho(Centers(i)),delta(Centers(i)),'o','MarkerSize',10,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
    hold on
end
if d==2
    if n<10000
        markerSize = 5;
    else
        markerSize = 1;
    end
    figure(2)
    for i=1:NCLUST
        ic=int8((i*64.)/(NCLUST*1.));
        plot(data((CL==i),1),data((CL==i),2),'o','MarkerSize',markerSize,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
        hold on
    end
    title ('Clustering result','FontSize',15.0);
end
end
