function [mrca, complete_genealogy, coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m)

generations = size(genealogy_m,1); %establish the dimentions of the genealogy_m matrix to set up the loops
lineages = 1:size(genealogy_m,2); %establish the dimentions of the genealogy_m matrix to set up the loops

coal_count = 0;
coal_events = -1*ones(length(lineages),3);

time_count = 0; %keeps track of the time, to establish t-mrca

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
            fecundities = round(age_dist_m(:,g-1).*transpose(leslie_matrix(1,:))); %vector of the expected number of offspring from parents of each age class
            parent_age_index = randi(sum(fecundities)); %randomly choose a number from the expected number of offspring
            for a = 1:length(fecundities) %SWITCH TO PROBABILITIES FROM EXPECTED VALUES?
                if parent_age_index <= sum(fecundities(1:a))
                    parent_age = a-1; %assign parent age based on randomly chosen parent index.
                    break
                end
            end
            
            genealogy_m(g-1,k,1) = randi([options_lower(parent_age+1) options_upper(parent_age+1)]); %assign a parent from the chosen age class
            genealogy_m(g-1,k,2) = parent_age; %set the age of the individual in the previous generation.
            
        else %case where the individual is not a newborn and age-1 ancestor must be chosen
            
            age_old = genealogy_m(g,k,2); %set an age variable equal to the current age in lineage k
            genealogy_m(g-1,k,1) = randi([options_lower(age_old) options_upper(age_old)]);
            genealogy_m(g-1,k,2) = age_old-1; %set the age of the individual in the previous generation.
            if k>1
                for r = 1:k-1
                    if (genealogy_m(g,r,1) ~= 0)
                    while (genealogy_m(g-1,r,2) == genealogy_m(g-1,k,2)) && (genealogy_m(g-1,r,1) == genealogy_m(g-1,k,1)) %consider the case where two individuals age back to the exact same individual in the previous generation
                        genealogy_m(g-1,k,1) = randi([options_lower(age_old) options_upper(age_old)]); %redraw a new individual if the exact same individual is drawn
                    end
                    end
                end
            end
        end
    end
    for r = 1:length(lineages-1)
        for s = r+1:length(lineages)
            if (isequal(genealogy_m(g-1,r,1), genealogy_m(g-1,s,1))) %if there is a coalescent event
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
end
complete_genealogy = genealogy_m;
mrca = time_count;

end