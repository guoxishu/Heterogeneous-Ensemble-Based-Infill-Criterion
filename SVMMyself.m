function models=SVMMyself(tr_x, tr_y)

global M
%tic;
[train_x,ps1]=mapminmax(tr_x',0,1);train_x=train_x';
[train_y,ps2]=mapminmax(tr_y',0,1);train_y=train_y';
models{1}=ps1;models{2}=ps2;
for i=1:M
    models{2+i} = svmtrain(train_y(:,i),train_x,'-s 3 -t 3 -q');
end
%toc;
%a


