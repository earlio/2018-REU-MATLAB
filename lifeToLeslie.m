function [LM] = lifeToLeslie(M)
%input should be an agex3 life table where col one is the age, two is the fecundity, two is the survival. 
%Output is corresponding Leslie Matrix

[x,y] = size(M);
LM = zeros(x,x);
for i = 1:x-1
   LM(1,i) = M(i,2);
   LM(i+1,i) = M(i,3);
end




