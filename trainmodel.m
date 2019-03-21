function models=trainmodel(tr_x, tr_y, index, str)

%PCA
[~, M]=PCAMyself(tr_x);
trx_fe=tr_x*M;

%FS
trx_fs=tr_x(:,index);

for i=1:length(str)
    str3=str(i,:);
    if strcmp(str3(1), 'FE')
        tr_xx=trx_fe; models(i).M=M;
    elseif strcmp(str3(1), 'FS')
        tr_xx=trx_fs; models(i).M=index;
    elseif strcmp(str3(1), 'NONE')
        tr_xx=tr_x; models(i).M=[];        
    end
    
    if strcmp(str3(2), 'RBF1')
        models(i).p=RBFMyself1(tr_xx, tr_y);
    elseif strcmp(str3(2), 'SVM')
        models(i).p=SVMMyself(tr_xx, tr_y);
    elseif strcmp(str3(2), 'RBF2')
        models(i).p=RBFMyself2(tr_xx, tr_y);
    end
end

        
        
        
        