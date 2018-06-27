function [indices] = terminal_indices(terminal_population, lineage_count, age_dist_m, age)
%returns a given number of indices in the terminal population

indices = zeros(2, lineage_count); %initializes the indices matrix which will fill the top surface of the genealogy matrix

%% Specified Same Age %%
if age > -1
for k = 1:lineage_count
    if isequal(age,0) %age 0 is a boundary case
        options = age_dist_m(end,1); %establish the possible indices to draw from, stored in a variable called options
        indices(1,k) = randi(options); %assigns an index from the possible individuals aged 0
        indices(2,:) = age; %assigns the second row of the indices matrix to 0
    else
        options = age_dist_m(end,age+1); %establishes the number of possible indices to draw from when age is not zero
        indices(1,k) = sum(age_dist_m(end,1:age)) + randi(options); %assigns an index from the possible individuals of specified age
        indices(2,:) = age; %assigns the second row of the indices matrix to the age
    end
end
        
%% Choose Two Random Individuals %%
else  
    
for k = 1:lineage_count
    indices(1,k) = randi(length(terminal_population)); %chooses random individuals from the terminal population
    indices(2,k) = terminal_population(indices(1,k)); %records the age corresponding to each index
end

end
