
%% Initial Conditions %%

pop_size = 50; %set the population size (fixed in the moran model)
number_generations = 50; %set the number of generations
lineage_count = 2; %set the number of lineages to trace

%% Choose Initial Lineages %%

indices = moran_indices(pop_size, lineage_count); %calls moran_indices to choose initial lineages to track

%% Track the Lineages %%

genealogy_moran = -1*ones(number_generations, lineage_count); %creates a genealogy matrix to track lineages
genealogy_moran(end,:) = indices; %sets the last row of the genealogy matrix equal to the initial indices

count = 0; %initializes generation counter
for i = number_generations:-1:1
    new_indices = choose_previous_indices(pop_size, lineage_count); %picks indices from the previous generation
    if isequal(new_indices(1),new_indices(2))
        genealogy_moran(i, 1) = new_indices(1); %coalesce to the first lineage if the 2 indices are the same
        break %if a coalescent event occurs, break the loop
    else
        genealogy_moran(i,:) = new_indices; %if there isn't a coalescent event, set the indices in the previous generation = new_indices
        count = count + 1; %add 1 to the generation counter
    end
end

disp(genealogy_moran); %displays the genealogy matrix
if isequal(count,50)
    disp('No MRCA'); %if no mrca is found, display no mrca
else
    mrca = count %if an mrca is found, display the count
end




