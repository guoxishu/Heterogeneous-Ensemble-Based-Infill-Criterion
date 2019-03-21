clear;clc;

itr=20;
global M V MaxValue MinValue 
Problems={'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ5','DTLZ6','DTLZ7',...
    'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'};
rand('seed', sum(100 * clock));
load ini10.mat

for Prob = 1:length(Problems)
    clear result
    Problem=Problems{Prob};
    k = find(~isstrprop(Problem,'digit'),1,'last');
    [V,M,iter,MaxValue,MinValue]=settings(Problem);
    for i=1:itr
        %% ≥ı ºevaluateµ„
        t1=clock;
        tr_x=ini(i).chrom;
        tr_x = tr_x.*repmat(MaxValue,11*V-1,1)+(1-tr_x).*repmat(MinValue,11*V-1,1);     
        tr_y=obj_real(tr_x,Problem);    
        [chrom, timestore, S]=nsga_opt(tr_x, tr_y, Problem, i, itr, iter);
        t=etime(clock, t1);
        result(i).ch=chrom;
        result(i).time=timestore;
        result(i).alltime=t;
        result(i).Scharacter=S;
    end
    a=[Problem(1) Problem(k+1) 'EN'];
    save(a,'result');
end