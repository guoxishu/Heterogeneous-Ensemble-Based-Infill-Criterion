function objv_LCB=estimate(x, tr_x, tr_y, models, str)

alpha=2;
TestSamNum=size(x,1);n=length(models);
sum=zeros(TestSamNum, size(tr_y,2));
for i=1:n
    clear y;
    stri=str(i,:);
    M=models(i).M;p=models(i).p;
    if strcmp(stri(1), 'FE')%PCA
        x_pca=x*M;
        if strcmp(stri(2), 'RBF1')%RBF
            y= sim(p,x_pca');y=y';
        elseif strcmp(stri(2), 'SVM')
            y=SVMtest(x_pca, p, TestSamNum);
        elseif strcmp(stri(2), 'RBF2')%RBF
            y=RBF2test(x_pca, p, TestSamNum);
        end
        
    elseif strcmp(stri(1), 'NONE')%NONE
        if strcmp(stri(2), 'RBF1')%RBF
            y= sim(p,x');y=y';
        elseif strcmp(stri(2), 'SVM')
            y=SVMtest(x, p, TestSamNum);
        elseif strcmp(stri(2), 'RBF2')%RBF
            y=RBF2test(x, p, TestSamNum);
        end
    
    elseif strcmp(stri(1), 'FS')%CSO
        x_cso=x(:,M);
        if strcmp(stri(2), 'RBF1')%RBF
            y= sim(p,x_cso');y=y';
        elseif strcmp(stri(2), 'SVM')
            y=SVMtest(x_cso, p, TestSamNum);
        elseif strcmp(stri(2), 'RBF2')%RBF
            y=RBF2test(x_cso, p, TestSamNum);
        end 
    end
    result(i).y=y;
    sum=sum+y;
end
me=sum/n;
sum=zeros(TestSamNum, size(tr_y,2));
for i=1:n
    y=result(i).y;
    sum=sum+(y-me).^2;
end
s2=sum/(n-1);
%fmin=min(tr_y);
%fmin_matrix=repmat(fmin,TestSamNum,1);
objv_LCB=me-alpha*sqrt(s2);

    


