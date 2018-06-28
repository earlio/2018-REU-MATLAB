function [indices] = terminal_indices(terminal_population, lineage_count)
%returns a given number of indices in the terminal population
indices = zeros(1,lineage_count);
for k = 1:lineage_count
    indices(k) = randi(length(terminal_population)); %chooses random individuals from the terminal population
end


