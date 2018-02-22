function [fit,err] = mymulti(x, ratio)
global array_basic_z;
global randomscore;
global h_score;

ind = find(x==1);
k = length(ind);

if (k == 0)               % empty module generated during deletion of node
    err = -1;
elseif (k == 1)           % single node module: algebraic connectivity would fail
    err = 0;
else
    err = connectivity_constrain(x);
end

if (k == 0)        
    f1 = 0;
    f2 = 0;
else
    % for Heinz score
    f1 = - sum(h_score(ind));
    
    % for normalized z score
%     f1 = - sum(array_basic_z(ind))/sqrt(k);    
    
    % for raw z score
%     f1 = -(sum(array_basic_z(ind))/sqrt(k) - randomscore(k,1))/randomscore(k,2);
%     f1 = -(sum(array_basic_z(ind))*randomscore(k,4) - randomscore(k,3));

    f2 = kegg_coverscore(x, ratio);
    
%     test one objectve
%     f2 = 0;
end

fit=[f1 f2];

