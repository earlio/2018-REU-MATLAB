function [mrca, complete_genealogy] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m)

generations = size(genealogy_m,1); %establish the dimentions of the genealogy_m matrix to set up the loops
lineages = size(genealogy_m,2); %establish the dimentions of the genealogy_m matrix to set up the loops

time_count = 0; %keeps track of the time, to establish t-mrca

for g = generations:-1:2 %iterate over the generations
    
    for k = lineages:-1:1 %iterate over the number of lineages
        
        if genealogy_m(g,k,2) == 0 %case where the individual is a newborn and a parent must be chosen
           fecundities = round(age_dist_m(:,g-1).*transpose(leslie_matrix(1,:)));
           parent_age_index = randi(sum(fecundities));
           for a = 1:length(fecundities)
               if parent_age_index <= sum(fecundities(1:a))
                   parent_age = a-1;
                   break
               end
           end
           choice = randi(age_dist_m(parent_age+1,g-1)); %establish choose a random number within the number of possible options for that age class
           genealogy_m(g-1,k,1) = sum(age_dist_m(1:parent_age,g-1)) + choice; %add the random number to the number of individuals in previous age classes to return an index for the previous generation
           genealogy_m(g-1,k,2) = parent_age; %set the age of the individual in the previous generation. 
           
        else %case where the individual is not a newborn and age-1 ancestor must be chosen
        %while (k==2) && (genealogy_m(g,k,2) ~= 0) && (genealogy_m(g-1,2,1) == genealogy_m(g-1,1,1))
           age_old = genealogy_m(g,k,2); %set an age variable equal to the current age in lineage k
           choice = randi(age_dist_m(age_old+1,g-1)); %establish choose a random number within the number of possible options for that age class
           if isequal(age_old, 1)
               genealogy_m(g-1,k,1) = choice;
               genealogy_m(g-1,k,2) = 0;
           else
           genealogy_m(g-1,k,1) = sum(age_dist_m(1:age_old,g-1)) + choice; %add the random number to the number of individuals in previous age classes to return an index for the previous generation
           genealogy_m(g-1,k,2) = age_old-1; %set the age of the individual in the previous generation. 
           end
        end
    end
time_count = time_count + 1; %add one to time count
if (isequal(genealogy_m(g-1,1,1), genealogy_m(g-1,2,1)))
    break
end
end

complete_genealogy = genealogy_m;
mrca = time_count;    