function fitness=PSO_CostFunction(p, x, y)

m=size(p,1);
fitness=100*ones(m,1);
for i=1:m
    index=find(p(i,:)==1);
    if ~isempty(index)
        xx=x(:,index);
        D=Discor(xx, y);
        R=Discor(xx,[]);
        fitness(i)=-0.8*D+0.2*R;
    end
end