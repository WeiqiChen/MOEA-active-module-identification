function child_offspring  = mutation_delete_node(parent_selected)
%% Description
% 1. Mutation by deleting node with lowest score
% 2. Input is in 'parent_selected' matrix of size [pop_size,V].
% 3. Output is also of same size in 'child_offspring'. 

%% Addition mutation
global V pop_size array_basic_z

T = max(array_basic_z) - min(array_basic_z);
child_offspring = zeros(pop_size,V);
for iteration = 1:pop_size
    mutated_child = parent_selected(iteration,:);  
    select_id = find(mutated_child);
   
    if (~isempty(select_id)) 
        select_zscore = array_basic_z(select_id);
        [minscore, index] = min(select_zscore);   % minimum score
        ind = select_id(index);          % the node with low z score
        prob = exp((min(array_basic_z) - minscore) / T);
        if (rand(1) < prob)            % stochastic
            mutated_child(ind) = 0;       %delete this node from mutated child
        end
    end 
   
    child_offspring(iteration,:) = mutated_child;
end