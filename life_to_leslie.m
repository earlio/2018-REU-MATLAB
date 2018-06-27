function [leslie_matrix] = life_to_leslie(file, sheet, cell_range)
%Input - a life table in a .csv or an .xlsx file in the form of an ages x 3 matrix. Column 1 contains the ages, Column 2 contains survival, Column 3 contains fecundity.  
%Output - corresponding Leslie Matrix
filename = file; %assign the file to a variable
if isequal(filename(end-2:end),'csv')
    life_table_m = csvread(filename,2,0); %read the csv file and assign it to a matrix called life_table_m. 2,0 offsets the entries to begin with teh values
elseif isequal(filename(end-2:end),'lsx')
    life_table_m = xlsread(filename, sheet, cell_range);
end

[ages] = size(life_table_m,1); %determine the number of rows of the life_table to create a correctly sized leslie matrix
leslie_matrix = zeros(ages,ages); %create a square leslie matrix with dimension = number of ages in the life table
for i = 1:ages-1 
   leslie_matrix(1,i) = life_table_m(i,2); %set the first row of the leslie matrix to the age-depdnent fecundity coefficients
   leslie_matrix(i+1,i) = life_table_m(i,3); %set the diagonal (beginning at (2,1)) equal to the age-dependent survival coefficients
end



