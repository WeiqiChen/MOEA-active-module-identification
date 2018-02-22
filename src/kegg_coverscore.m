function coverscore_minus = kegg_coverscore(x, ratio)
global KEGG_node KEGG_pos V 
global KEGG_node_colsize;

% k = length(index_subnet_monte);
% index_subnet_monte = find(x == 1);

%% maximun coverage for one single pathway
% maxp = 0;
% % maxpid = 0;
% for i = 1:size(KEGG_node,2)
%     covrate = length(intersect(KEGG_pos{i},index_subnet_monte));
%     if covrate > maxp
%         maxp = covrate;
% %         maxpid = i;
%     end
% end
% coverscore_minus = -covrate;
% % tmp(1) = maxp;
% % tmp(2) = maxpid;

%% total coverage throughout all pathways
% sum = 0;
% for i = 1:size(KEGG_node, 2)
%     sum = sum + length(intersect(KEGG_pos{i},index_subnet_monte));
% end
% 
% coverscore_minus = -sum/k;


%% sum of pathway cover rate through all pathways, vectorized version
count = x * KEGG_node;
threshold = ratio * KEGG_node_colsize;
coverscore_minus = - sum(count >= threshold);

%% unvectorized
% ind = find(x==1);
% count = zeros(1, size(KEGG_node, 2));
% for i = 1:size(KEGG_node, 2)
%     count(i)= length(intersect(KEGG_pos{i},ind));
% end
% threshold = ratio * KEGG_node_colsize;
% coverscore_minus = - sum(count >= threshold);

 %% test time
% tic
% for i = 1:1000
% %     count = x(i,:) * KEGG_node;
%     
%     ind = find(x(i,:)==1);
%     count = zeros(1, size(KEGG_node, 2));
%     for j = 1:size(KEGG_node, 2)
%         count(j)= length(intersect(KEGG_pos{j},ind));
%     end
% end
% toc


    