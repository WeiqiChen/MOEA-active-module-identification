function population = clearing(population)

global V M 
psize = size(population, 1);
population = sortrows(population, V+1);
markrow = population(1, 1:V);
for i = 2:psize
    currentrow = population(i, 1:V);
    if isequal(markrow, currentrow)  % find replicate of mark row
        population(i, V+M+1) = -2;
    else
        markrow = currentrow;
    end
end