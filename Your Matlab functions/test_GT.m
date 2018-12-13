function test_GT()

close all
clc;

   options = struct();
tic
[ETA DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = GT(1e5,options,0);
toc
if length(ETA) == 6
   fprintf('\n eta have the correct size') 
else
   fprintf('\n eta have not the correct size')    
end

try
    test = ETA(round(1+5*rand()));
    if test<    =1&&test>=0
        fprintf('\n efficiency correct order of magnitude') 
    else
        fprintf('\n efficiency higher than 1 or lower than 0')  
    end
catch
     fprintf('\n Problem with efficiencies...')  
end


if length(DATEN) == 2
   fprintf('\n DATEN have the correct size') 
else
   fprintf('\n DATEN have not the correct size')    
end

if length(DATEX) == 4
   fprintf('\n DATEX have the correct size') 
else
   fprintf('\n DATEX have not the correct size')    
end

if (size(DAT,1) == 5) && (size(DAT,2) == 4)
   fprintf('\n DAT have the correct size') 
else
   fprintf('\n DAT have NOT the correct size')    
end


try
   isfield(COMBUSTION,'LHV') ;
   isfield(COMBUSTION,'e_c') ;
   isfield(COMBUSTION,'lambda') ;
   isfield(COMBUSTION,'Cp_g') ;
   isfield(COMBUSTION,'fum') ;
   fprintf('\n Combustion have the right structure');
catch
   fprintf('\n Combustion have not the right structure');    
end


end