function [indices] = terminal_indices(lineage_count, age_dist_m, age)
%returns a given number of indices in the terminal population

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
    %indices(1,k) = options_lower(age+1) + randi(age_dist_m(age+1,generations)); %set the index
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
                indices(1,k) = options_lower(age+1) + round((options_upper(age+1)-options_lower(age+1))*rand); %reset the index if it has already been used
                indices(2,k) = age; %set the age
            end
        end
    end
    age = stored_age;
end

