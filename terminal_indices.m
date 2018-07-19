function [indices] = terminal_indices(lineage_count, age_dist_m, age)
<<<<<<< HEAD
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
=======
%Inputs:
%1. lineage_count - the preset number of lineages that the user wants to
%sample
%2. age_dist_m - the age distribution/demographic matrix containing the
%population of each age class over the number of generations
%3. age - should be -1 if #lineage_count random individuals are to be
%selected, or a given age >= 0 if #lineage_count individuals from that
%specific age class are to be sampled.

%This function randomly selects indices in the final generation of the
%demographic matrix to begin the lineage trace. If the age input is -1 it
%randomly selects #lineage_count individuals randomly (without
%replacement). If the age input is 0 or higher it randomly selects
%individuals from that given age class

%Outputs:
%1. indices - indices is a 2xlineage_count matrix which contains both the
%index of each selected individual and its age. The first row contains the
%indices and the second row contains the ages.

indices = zeros(2, lineage_count); %initializes the indices matrix which will fill the top surface of the genealogy matrix
[age_classes,generations] = size(age_dist_m); %measures size of age_dist_m
options_lower = [1]; %sets a lower bound vector for each age class
options_upper = [age_dist_m(1,end)]; %sets an upper bound vector for each age class
for a = 2:size(age_dist_m)
    options_lower = [options_lower sum(age_dist_m(1:a-1,end))+1]; %construct a vector of lower bounds
    options_upper = [options_upper sum(age_dist_m(1:a,end))]; %construct a vector of upper bounds
end
scaling = sum(age_dist_m(1:age_classes,generations)); %scale factor for age class sizes
probabilities = age_dist_m(1:age_classes,generations)./scaling; %scales the age class sizes as probabilities
stored_age = age; %stores the age for each iteration of the loop

%% Specified Same Age %%
for k = 1:lineage_count
    if isequal(age,-1)
        age_index = rand; %chooses random number between zero and one
        for a = 1:age_classes
            if age_index <= sum(probabilities(1:a)) %check if the random number is less than the values of the previous probabilities combined
                age = a-1; %if the random number falls in the corresponding range of an age class, assign the individual to that age class
                break %the loop is broken when an age is assigned
            end
        end
    end
    indices(1,k) = options_lower(age+1) + round((options_upper(age+1)-options_lower(age+1))*rand);
    indices(2,k) = age;
    if k>1
        for r = 1:k-1
            while isequal(indices(1,k),indices(1,r))
                age_index = rand; %chooses random number between zero and one
                for a = 1:age_classes
                    if age_index <= sum(probabilities(1:a)) %check if the random number is less than the values of the previous probabilities combined
                        age = a-1; %if the random number falls in the corresponding range of an age class, assign the individual to that age class
                        break %the loop is broken when an age is assigned
                    end
                end
                indices(1,k) = options_lower(age+1) + round((options_upper(age+1)-options_lower(age+1))*rand);
                indices(2,k) = age;
                
            end
        end
    end
    age = stored_age;
>>>>>>> L2Lv5
end

