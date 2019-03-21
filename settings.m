function [V,M,gen,MaxValue,MinValue]=settings(Problem)

    V = 10;gen=24;
    switch Problem
        case {'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7'}
            M = 3;MaxValue = ones(1,V);MinValue = zeros(1,V);
            
        case {'ZDT6'}
            M = 2;MaxValue = ones(1,V);MinValue = zeros(1,V);   
            
        case {'ZDT3'}   
            M = 2;MaxValue = ones(1,V);MinValue = zeros(1,V);
            
        case {'ZDT4'}
            M = 2;MaxValue = [ones(1,M-1) 5*ones(1,V-M+1)];MinValue = [zeros(1,M-1) -5*ones(1,V-M+1)];            
        
        case {'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'} 
            M = 3;MaxValue = 2*(1:V);MinValue = zeros(1,V);           
    end

end
