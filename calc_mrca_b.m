function [mrca, complete_genealogy, coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m)

%Inputs:
%1. genealogy_m - a number_generations x lineage_count x 2 matrix. At the
%time of input, the matrix should be filled with a null value of -1 at
%every index except it should be filled with the indices and ages of
%individuals in it's final surface/row.
%2. leslie_matrix -  a leslie matrix corresponding to the input life table,
%generated using the life_to_leslie function
%3. age_dist_m - the age distribution/demographic matrix containing the
%population of each age class over the number of generations

%this function traces a preset number of lineages (lineage_count) back in
%time to find a most recent common ancestor. The function is written
%iteratively and it's loops are broken if a common ancestor is reached. 

%Outputs: 
%1. mrca - the time to most recent common ancestor. Mrca is equal to the maximum
%number of generations if no most recent common ancestor is reached. 

%2. complete_genealogy - a completed genealogy_m. Page one contains the indices of
%individuals in each lineage for every generation until a most recent
%common ancestor is reached. Page two contains the ages of individuals in
%each lineage for every generation until a most recent common ancestor is
%reached. When an MRCA is reached, the lineages will converge and all
%remaining cells will be filled with "-1".


%3. coal_events - coal events is a vector (when lineage_count = 2) or a matrix (when lineage_count > 2)  
% where the first column is records the lineage lost in a coalescent event,
% the second column is the remaining lineage after a coalescent event, and
% the third column is the generation count when the coalescent event
% occurs. 

generations = size(genealogy_m,1); %establish the dimentions of the genealogy_m matrix to set up the loops
lineages = 1:size(genealogy_m,2); %establish the dimentions of the genealogy_m matrix to set up the loops

coal_count = 0;
coal_events = -1*ones(length(lineages)-1,3);

time_count = 1; %keeps track of the time, to establish t-mrca

for g = generations:-1:2 %iterate over the generations
    options_lower = [1]; %sets a lower bound vector for each age class
    options_upper = [age_dist_m(1,g-1)]; %sets an upper bound vector for each age class
    for a = 2:size(age_dist_m)
        options_lower = [options_lower sum(age_dist_m(1:a-1,g-1))+1]; %construct a vector of lower bounds
        options_upper = [options_upper sum(age_dist_m(1:a,g-1))]; %construct a vector of upper bounds
    end
    for q = 1:length(lineages) %iterate over the number of lineages
        k = lineages(q);
        if genealogy_m(g,k,2) == 0 %case where the individual is a newborn and a parent must be chosen
            fecundities = age_dist_m(:,g-1).*transpose(leslie_matrix(1,:)); %great a vector of expected number of offspring
            scaling = sum(fecundities); % determines the value to scale the fecundities by
            fecundities = fecundities./scaling; %scales the fecundities vector so the sum is one
            parent_age_index = rand; %chooses a uniformly dist random number between 0 and 1
            for a = 1:length(fecundities)
                if parent_age_index <= sum(fecundities(1:a)) %check if the random number is less than the values of the previous fecundities combined
                    parent_age = a-1; %if so the parent age is assigned to that age value
                    break %the loop is broken when a parent age is assigned
                end
            end
%             if isequal(age_dist_m(parent_age+1,g-1),1) %check to ensure there is more than one individual in parent age class
%             for p = k:length(lineages)
%                 j = lineages(p); %check all remaining lineages > k
%                 if isequal(parent_age, (genealogy_m(g,j,2)-1)) 
%                     parent_age = parent_age -1; %decrease parent age if entering an age class with one individual but an individual from a different lineage is going to age into that age class
%                     break
%                 end
%             end
%             end
            genealogy_m(g-1,k,1) = options_lower(parent_age+1) + round((options_upper(parent_age+1)-options_lower(parent_age+1))*rand); %assign a parent from the chosen age class
            genealogy_m(g-1,k,2) = parent_age; %set the age of the individual in the previous generation.
            
        else %case where the individual is not a newborn and age-1 ancestor must be chosen
            
            age_old = genealogy_m(g,k,2); %set an age variable equal to the current age in lineage k
            genealogy_m(g-1,k,2) = age_old-1; %set the age of the individual in the previous generation.
            genealogy_m(g-1,k,1) = options_lower(age_old+1) + round((options_upper(age_old+1)-options_lower(age_old+1))*rand);
            if k>1    
                for r = 1:k-1
                    if isequal(genealogy_m(g-1,r,2),0)
                        break
                    else
                    while isequal(genealogy_m(g,r,2),genealogy_m(g,k,2)) && isequal(genealogy_m(g-1,r,1),genealogy_m(g-1,k,1))
                        genealogy_m(g-1,k,1) = options_lower(age_old+1) + round((options_upper(age_old+1)-options_lower(age_old+1))*rand);
                    end               
                    end
                end
            end
        end
    end
    for r = 1:length(lineages-1)
        for s = r+1:length(lineages)
            if (isequal(genealogy_m(g-1,r,1), genealogy_m(g-1,s,1))) %if there is a coalescent event
               % genealogy_m(g-1,s,:) = -1;
                lineages(s)=[]; %remove the later lineage
                coal_count = coal_count+1;
                coal_events(coal_count,:) = [s r time_count];
            end
        end
    end
time_count = time_count + 1; %add one to time
if isequal(length(lineages),1)
    break %breaks the loop if the number of lineages reaches 1(at an mrca)
end
if isequal(time_count,generations)
    break
end
end
complete_genealogy = genealogy_m;
mrca = time_count;