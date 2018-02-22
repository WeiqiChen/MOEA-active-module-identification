
%% Code Starts
%% read in initial files
% read in p value
% array_p_value = textread('CDpvalue.txt','%f');
array_p_value = textread('pvalue.txt', '%f');

% read in triple matrix and convert to adjacency matrix
% tri_matrix = textread('CDtriplematrix.txt');
tri_matrix = textread('triplematrix.txt');
n = length(array_p_value);
G = sparse([tri_matrix(:,1);tri_matrix(:,2)], [tri_matrix(:,2);tri_matrix(:,1)],...
    [tri_matrix(:,3), tri_matrix(:,3)], n, n);

% read in KEGG map
% KEGG_node = textread('CD_KEGG_map.txt');
KEGG_node = textread('kegg_map.txt');

%% JActiveModule method
%% calculate z score from p value
array_basic_z = getZscore(array_p_value);

% %% Normalization
% mu = mean(array_basic_z);
% sigma = std(array_basic_z);
% array_basic_z = (array_basic_z - mu) ./ sigma;

%% Test z score 
z_sort = sort(array_basic_z,'descend');
num = length(z_sort);
lambda = [0, 0.001, 0.003, 0.01, 0.03, 0.1];
mscore = zeros(1, num);
for j = 1:length(lambda)
    for i = 1:num
%         normalized z score
%         mscore(i) = sum(z_sort(1:i))/sqrt(i) - lambda(j)*i;
%         mscore(i) = (sum(z_sort(1:i))/sqrt(i) - randomscore(i, 1)) / randomscore(i, 2) - lambda(j)*i;
%         un_normalized z score
        mscore(i) = (sum(z_sort(1:i))/sqrt(i) - randomscore(i, 1)) / randomscore(i, 2) - lambda(j)*i;
%         mscore(i) = sum(z_sort(1:i))/sqrt(i) - lambda(j)*i;
    end   
    subplot(2, 3, j);
    plot(mscore);
    type = ['background correction, lambda = ',num2str(lambda(j))];
    title({'module score on un-normalized data (original)';type});
end


%% perform monte carlo approach

iter_monte = 10000;
graphSize = length(G);
randomscore = zeros(graphSize, 2);
for i = 1:graphSize
    [mu,sigma] = monte(graphSize,i,array_basic_z,iter_monte);
    randomscore(i,1) =  mu;
    randomscore(i,2) =  sigma;
    disp(['monte carlo simulation: completed ' num2str(i) ' of ' num2str(graphSize) ' iterations']);
end

% Pre-compute: column 3 as mu/sigma, column 4 as 1/(sigma*sqrt(k))
mu2sigma = randomscore(:,1)./randomscore(:,2);
column4 = zeros(graphSize, 1);
for k = 1:graphSize
    column4(k) = 1/(randomscore(k, 2)*sqrt(k));
end

randomscore = horzcat(randomscore, mu2sigma, column4);

%% Heinz method
%% Test Heinz method 
tau = [0.000003, 0.00001, 0.00003, 0.0001, 0.0003, 0.001];
sorted_p = sort(array_p_value);
num = length(array_p_value);
mscore = zeros(1, num);
    
for j = 1:length(tau)
    h_score = log(tau(j)) - log(sorted_p);   
    for i = 1:num
       mscore(i) = sum(h_score(1:i));
    end
    subplot(2, 3, j);
    plot(mscore);
    type = ['threshold = ',num2str(tau(j))];
    title({'Heinz module score on p value';type});
end

%% Estimate parameters for beta distribution and test
para = mle(array_p_value,'distribution','beta');
alpha = para(1);
beta = para(2);
y = betarnd(alpha, beta, 330, 1);
subplot(1,3,1)
hist(array_p_value, 50)
subplot(1,3,2)
hist(y, 50)
subplot(1,3,3)
qqplot(array_p_value, y)

%% Test parameters estimated from BUM model using BioNet R package

lambda = 0.09066589;
a = 0.1133399;

% tau is calculated under the FDR alpha = 0.001; maybe changed
tau = 0.002366535;
% x = 0.001:0.001:1;
% y = lambda * rand(1, 1000) + (1-lambda)*a.*(x.^(a-1));
% y = lambda * rand(330, 1) + (1-lambda)* betarnd(a, 1, 330, 1);
y = betarnd(a, 1, 330, 1);
subplot(1,3,1)
hist(array_p_value, 50)
subplot(1,3,2)
hist(y, 50)
subplot(1,3,3)
qqplot(array_p_value, y)


%% processing KEGG data

KEGG_node_colsize = sum(KEGG_node,1);   % number of genes in each pathway

KEGG_pos = cell(1, size(KEGG_node, 2));
for i = 1:length(KEGG_pos)
    KEGG_pos{i} = find(KEGG_node(:,i));
end

KEGG_gene = sum(KEGG_node, 2);    % number of pathways each gene is involved
%% save as mat data file
save yeast_drug.mat G array_p_value KEGG_node KEGG_pos KEGG_node_colsize KEGG_gene a tau
clear

%% Pre-processing ends; ready for community detection


%% Merge figures

% Load saved figures
% c=hgload('MyFirstFigure.fig');
% k=hgload('MySecondFigure.fig');
list = [3, 5, 7, 9];
out = zeros(1,4);
for ii = 1:length(list)
    i = list(ii);
    in = strcat('300pop, 2000gen, 0.0005tau, 0.',num2str(i),'ratio.fig');
    out(i) = hgload(in);
end

% Prepare subplots
figure
% h(1)=subplot(1,2,1);
% h(2)=subplot(1,2,2);
for ii = 1:length(list)
    i = list(ii);
    h(i) = subplot(3, 3,i);
end

% Paste figures on the subplots
% copyobj(allchild(get(c,'CurrentAxes')),h(1));
% copyobj(allchild(get(k,'CurrentAxes')),h(2));
for ii = 1:length(list)
    i = list(ii);
    copyobj(allchild(get(out(i),'CurrentAxes')),h(i));
end

% Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')
% for ii = 1:length(list)
%     i = list(ii);
%     l(i)=legend(h(i),'LegendForFigure')
% end
