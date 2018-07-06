function[new_indices] = choose_previous_indices(pop_size, lineage_count)
%input is a fixed population size and a number of lineages to trace
%output is two indices corresponding to individuals in the previous
%generation

options = 1:pop_size; %create a vector of indices
removed_individual = randi(pop_size); %choose one individual to remove
parent = randi(pop_size); %choose one individual to be a parent
options(removed_individual) = parent; %replace the removed individual's index with the parent's index
new_indices = zeros(1,lineage_count); %create a vector to store the new indices. 
for k = 1:lineage_count
    choice = randi(length(options)); %pick a random index of the options vector
    new_indices(k) = options(choice); %choose random numbers and find the corresponding index for each lineage
    options(choice) = []; %remove that number from the options vector
end
    