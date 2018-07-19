function [indices] = moran_indices(pop_size, lineage_count)
%inputs are the moran population size and the number of lineages
%outputs are initial indices

indices = zeros(1,lineage_count); %establish an indices vector
options = 1:pop_size; %create a vector of index options
for i = 1:lineage_count
    choice = randi(length(options)); %choose a location on the options vector
    indices(i) = options(choice); %return that index to the indices vector
    options(choice) = []; %remove that index from the options vector
end

    