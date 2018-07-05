function [indices] = terminal_indices(lineage_count, age_dist_m, age)
%returns a given number of indices in the terminal population

indices = zeros(2, lineage_count); %initializes the indices matrix which will fill the top surface of the genealogy matrix
age_classes = size(age_dist_m,1);
%% Specified Same Age %%
if age > -1
for k = 1:lineage_count
    if isequal(age,0) %age 0 is a boundary case
        options = 1:age_dist_m(end,1);
        option_1 = randi(length(options));
        indices(1,k) = option_1; %assigns an index from the possible individuals aged 0
        options(option_1) = [];
        indices(2,:) = age; %assigns the second row of the indices matrix to 0
    else
        options = sum(age_dist_m(end,1:age)):sum(age_dist_m(end,1:age+1)); %establishes the number of possible indices to draw from when age is not zero
        option_2 = randi(length(options));
        indices(1,k) = option_2; %assigns an index from the possible individuals of specified age
        options(option_2) = [];
        indices(2,:) = age; %assigns the second row of the indices matrix to the age
    end
end
        
%% Choose Two Random Individuals %%
else  

choices = 1:sum(age_dist_m(:,end));

for k = 1:lineage_count
    slot_n = randi(length(choices));
    indices(1,k) = choices(slot_n); %chooses random individuals from the terminal population
    choices(slot_n) = [];
    for a = 1:age_classes
        if (sum(age_dist_m(1:a-1,end)) <= indices(1,k)) && (indices(1,k)<= sum(age_dist_m(1:a,end))) %records the age corresponding to each index
        indices(2,k) = a-1;
        break
        end
    end
end
end
