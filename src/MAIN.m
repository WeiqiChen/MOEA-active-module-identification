%% Description

% 1. This is the main program of NSGA II.
% 2. This code defines population size in 'pop_size', number of design
%    variables in 'V', number of runs in 'no_runs', maximum number of 
%    generations in 'gen_max', current generation in 'gen_count' and number of objectives
%    in 'M'.
% 3. Final optimal Pareto soutions are in the variable 'pareto_rank1', with design
%    variables in the coumns (1:V), objectives in the columns (V+1 to V+M),
%    constraint violation in the column (V+M+1), Rank in (V+M+2), Distance
%    in (V+M+3).




%% Code starts
%% load data
clear
clc
global V M pop_size ;

global G;
% global array_basic_z;
% global randomscore;
global KEGG_node;
global KEGG_pos;
global KEGG_node_colsize;
global KEGG_gene;

global h_score;

% load('yeastData.mat');
load('yeast_drug.mat');
G = full(G);

%% Parameters Initialization
M=2;                   % Two objectives
V=size(G,1);           % Variable number is the number of nodes in the whole network
pop_size=100;           % Population size
no_runs=1;              % Number of runs
gen_max=20;            % MAx number of generations - stopping criteria
fname='mymulti';        % Objective function and constraint evaluation
% Deletion_Probability = 0.3;

% defualt tau is calculated under the FDR alpha = 0.001
% set tau a different value (usually smaller)using following line
% tau = 1.5152e-007;   % for Heinz score calculation
h_score = (a - 1) * (log(array_p_value) - log(tau));  

ratio = [0.9];  % ratio for counting kegg pathway numbers
%% run starts
Q=[];
for run = 1:length(ratio)
% for run = 1:no_runs
    
   
%% Initial population 

% x = rand(pop_size,V)>0.5;
% [nComponents,sizes,members] = networkComponents(G);
% valid_node = members{1};
% choose = rand(pop_size,length(valid_node))>0.5;
% x = zeros(pop_size,V);
% for i = 1:pop_size
%     id = find(choose(i,:));
%     x(i,valid_node(id)) = 1;
% end

[sortedValues,indexValue] = sort(h_score,'descend');  % Sort the values in descending order
[sortedDeg,indexDeg] = sort(sum(G, 2),'descend'); 

sortIndex = [indexValue(1:(0.6*pop_size)); indexDeg(1:(0.4*pop_size))];
x = zeros(pop_size,V);

% Initialize population from one node
% If the first selected node is isolated, then only one node is initialized
k = length(h_score);
for i = 1:pop_size
    firstID = mod(i, k);
    if (firstID  == 0)
        firstID = V;
    end
    firstID = sortIndex(firstID);
    x(i,firstID) = 1;
    neighborID = find(G(firstID,:));
%     if (~isempty(neighborID))
%         secondID = neighborID(randi(length(neighborID)));
%         x(i,secondID) = 1;
%     end   
    ind = neighborID(h_score(neighborID) > 0);   % add all neighbor id with positive h score
    x(i, ind) = 1;
end


%% Evaluate objective function
ff = zeros(pop_size,M);
err = zeros(pop_size,1);
for i =1:pop_size
    [ff(i,:),err(i,:)] =feval(fname, x(i,:), ratio(run));           % Objective function evaulation 
%     if err(i,:) <= 0
%         x(i,:) = repair_disconnected(x(i,:));           % Repair infeasible solution
%         [ff(i,:),err(i,:)] =feval(fname, x(i,:));       % Get new objective functions for repaired solution
%     end
end
population_init=[x ff err];
[population,front]=NDS_CD_cons(population_init);    % Non domination Sorting on initial population
    
%% Generation Starts
for gen_count=1:gen_max
%% selection (Parent Pt of 'N' pop size)
% parent_selected=tour_selection(population);                     % 10 Tournament selection
parent_selected=population;    % Reduce the selection pressure to avoid premature

%% Reproduction (Offspring Qt of 'N' pop size)
% child_offspring  = genetic_operator(parent_selected(:,1:V));    % SBX crossover and polynomial mutation
% child_offspring  = cross_over(parent_selected(:,1:V));  %single point crossover and bitwise mutation


% child_offspring  = parent_selected(:,1:V);
child_offspring  = cross_over(parent_selected(:,1:V)); 

child_offspring  = mutation_add_node(child_offspring);
% if mod(gen_count, 2) == 1
%     child_offspring  = mutation_add_node(child_offspring);
% else
%     child_offspring  = mutation_delete_node(child_offspring);
% end
%%

for ii = 1:pop_size
    try
        [ff(ii,:),err(ii,:)] = feval(fname, child_offspring(ii,:), ratio(run));  % objective function evaluation for offspring
%         if err(ii,:) <= 0
%             child_offspring(ii,:) = repair_disconnected(child_offspring(ii,:));
%             [ff(ii,:),err(ii,:)]=feval(fname, child_offspring(ii,:));
%         end
    catch
        [ff(ii,:),err(ii,:)]=mymulti_fail_to_converge(child_offspring(ii,:), ratio(run));
    end
end

% error_norm=normalisation(err);                                  
% child_offspring=[child_offspring fff error_norm];
child_offspring=[child_offspring ff err];

%% INtermediate population (Rt= Pt U Qt of 2N size)
population_inter=[population(:,1:V+M+1) ; child_offspring(:,1:V+M+1)];

% a clearing procedure as a niching method
population_inter = clearing(population_inter);  

[population_inter_sorted,front]=NDS_CD_cons(population_inter);              % Non domination Sorting on offspring
%% Replacement - N
new_pop=replacement(population_inter_sorted, front);
population=new_pop;
disp(['run ' num2str(run) ': ' num2str(gen_count) ' in ' num2str(gen_max) ' generations']);
if mod(gen_count, 5) == 0
%     plot(new_pop(:,V+1),new_pop(:,V+2))
    plot(new_pop(:,V+1),new_pop(:,V+2),'*')
    hold on
end

% plot(new_pop(:,V+1),new_pop(:,V+2))
% plot(new_pop(:,V+1),new_pop(:,V+2),'*')
% hold on


end
new_pop=sortrows(new_pop,[V+M+2, V+1]);
paretoset(run).trial=new_pop(:,1:V+M+1);
Q = [Q; paretoset(run).trial];                      % Combining Pareto solutions obtained in each run


filename = ['NSGA2_Heinz_',num2str(pop_size),'pop_',num2str(gen_max), 'gen_',num2str(tau),...
    'tau_',num2str(ratio(run)),'ratio_CROSSOVER.csv'];
% filename = ['NSGA2_Heinz_',num2str(pop_size),'pop_',num2str(gen_max), 'gen_',num2str(tau),...
%     'tau_',num2str(ratio(run)),'ratio.csv'];
csvwrite(filename,new_pop);

figname = [num2str(pop_size),'pop, ',num2str(gen_max), 'gen, ',num2str(tau),...
    'tau,',num2str(ratio(run)),'ratio,CROSSOVER'];
% figname = [num2str(pop_size),'pop, ',num2str(gen_max), 'gen, ',num2str(tau),...
%     'tau,',num2str(ratio(run)),'ratio'];
title({'Evolved Fitness per 5 Generation'; figname});
xlabel('Module Score')
ylabel('KEGG Coverage Score')

saveas(1, [figname,'.fig'])
close all

end

%% Result and Pareto plot
if run==1
plot(new_pop(:,V+1),new_pop(:,V+2),'*')
else                                        
[pareto_filter,front]=NDS_CD_cons(Q);               % Applying non domination sorting on the combined Pareto solution set
rank1_index=find(pareto_filter(:,V+M+2)==1);        % Filtering the best solutions of rank 1 Pareto
pareto_rank1=pareto_filter(rank1_index,1:V+M);
plot(pareto_rank1(:,V+1),pareto_rank1(:,V+2),'*')   % Final Pareto plot
end

% title('Evolved Fitness per 200 Generation');
% xlabel('Active Module Score')
% ylabel('KEGG Coverage Score')
% str = ['One objective Test, ',num2str(gen_max),' generations'];
% str = [num2str(gen_max),' generations'];
% title(str );

% % disp('row and col sum of new population');
% (sum(new_pop(:,1:V),2))'
% sum(new_pop(:,1:V),1)

% filename = strcat('NSGA2_',num2str(pop_size),'pop_',num2str(gen_max),'gen_',num2str(Deletion_Probability),'delete_NonTour_RandomGrow.csv');
% filename = strcat('NSGA2_Heinz_',num2str(pop_size),'pop_',num2str(gen_max),...
%         'gen_',num2str(tau),'tau_NonTour_RandomGrow_SADelete.csv');
% csvwrite(filename,new_pop)

% colsum = sum(new_pop(:,1:V),1);
% favid = find(colsum);
% [sortzscore,sortzid] = sort(array_basic_z,'descend');
% topzid = sortzid(1:size(favid,2));
% overlap_id = intersect(favid, topzid);
% disp('array_basic_z of module id');
% array_basic_z(favid)'
% disp('array_basic_z of top z score');
% array_basic_z(topzid)'



