function child_offspring  = cross_over(parent_selected)
%% Description
% 1. Crossover followed by mutation
% 2. Input is in 'parent_selected' matrix of size [pop_size,V].
% 3. Output is also of same size in 'child_offspring'. 

%% 
global pop_size V
rownum = 0.5*pop_size;
rc = randi(pop_size, [rownum, 2]);
child_offspring = parent_selected;
%% no crossover, only mutation
% child_offspring = zeros(pop_size,V);
% for i=1:pop_size
%     parent = parent_selected(i,:);  
%     child_offspring(i,:) = bit_mutation(parent);           % bitwise mutation
% end

%% Single-point Cross over         % would generate disconnected modules
for i=1:size(rc, 1)
    parent1=parent_selected(rc(i, 1),:);  %Randomly choose 2 individuals
    parent2=parent_selected(rc(i, 2),:);
    cp = randi(V);                    %Choose a crossover point
    child1 = [parent1(1:cp), parent2(cp+1:end)]; %Crossover
    child2 = [parent2(1:cp), parent1(cp+1:end)];
%     child_offspring((rc(2*i-1)),:)=bit_mutation(child1);           % bitwise mutation
%     child_offspring((rc(2*i)),:)=bit_mutation(child2); 
    child_offspring(rc(i, 1),:) = child1;           
    child_offspring(rc(i, 2),:) = child2; 
end