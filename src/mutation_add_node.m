function child_offspring  = mutation_add_node(parent_selected)
%% Description
% 1. Crossover followed by mutation
% 2. Input is in 'parent_selected' matrix of size [pop_size,V].
% 3. Output is also of same size in 'child_offspring'. 

%% Addition mutation
global G V pop_size h_score
global KEGG_gene
% global array_basic_z

child_offspring = zeros(pop_size,V);
for iteration = 1:pop_size
    mutated_child = parent_selected(iteration,:);   
    select_id = find(sum(G(find(mutated_child),:)));  % get all the neighbor id
   
% %     select_zscore = array_basic_z(select_id);
% %     [sortValue,sortIndex] = sort(select_zscore,'descend'); %sort by descending order
%     for i=1:length(select_id)
% %        ind = select_id(sortIndex(i));  % select node by z score ranking from high to low
%        ind = select_id(randi(length(select_id)));   % select node randomly, avoid local trap, results in slow evolution
%        if mutated_child(ind) == 0
%            mutated_child(ind) = 1;        %bit mutation    
%            break
%        end
%     end
   
%     candidates = setdiff(select_id, parent_id);  % get neighboring nodes that are not in parent id
%     if (~isempty(candidates))
%         ind = candidates(randi(length(candidates)));  % randomly choose a node to add
%         mutated_child(ind) = 1;
%     end
    
    
%     if rand < 0.5
%         mutated_child(select_id) = 1;    % add all neighboring nodes
%     else
%         ind = select_id(h_score(select_id) > 0);   % add all neighboring nodes with positive score
%         mutated_child(ind) = 1;
%         ind = select_id(KEGG_gene(select_id) > 0);   % add all neighboring nodes in pathways
%         mutated_child(ind) = 1;
%     end
    
    ind = select_id(h_score(select_id) > 0);   % add all neighboring nodes with positive score
    mutated_child(ind) = 1;
    if rand < 0.5
        ind = select_id(KEGG_gene(select_id) > 1);   % add neighboring nodes in pathways
        mutated_child(ind) = 1;
    end
    child_offspring(iteration,:) = mutated_child;
    
end