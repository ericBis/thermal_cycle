function test_ST()

close all
clc;

   options = struct();
tic
[ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(1e5,options,0);
toc
if length(ETA) == 9
   fprintf('\n eta have the correct size') 
else
   fprintf('\n eta have not the correct size')    
end

try
    test = ETA(round(1+8*rand()));
    if test<=1&&test>=0
        fprintf('\n efficiency correct order of magnitude') 
    else
        fprintf('\n efficiency higher than 1 or lower than 0')  
    end
catch
     fprintf('\n Problem with efficiencies...')  
end

try
    for i=1:length(XMASSFLOW)
        if XMASSFLOW(i)<=0
            fprintf('\n Negative massflow :-/...');
        end
    end
catch
     fprintf('\n No feedheaters? Can be ok, if it is not imposed in the setup')  
end


if length(DATEN) == 3
   fprintf('\n DATEN have the correct size') 
else
   fprintf('\n DATEN have not the correct size')    
end

if length(DATEX) == 7
   fprintf('\n DATEX have the correct size') 
else
   fprintf('\n DATEX have not the correct size')    
end

fprintf('\n the output matrix DAT have a size of (%d,%d)',size(DAT,1),size(DAT,2));


if length(MASSFLOW) == 4
   fprintf('\n DATEN have the correct size') 
else
   fprintf('\n DATEN have not the correct size')    
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