function [chrom, timestore, selection]=nsga_opt(tr_x, tr_y, Problem, ceng, itr, iter)
    
    global M V %L1
    L=11*V-1+25;
    timestore=zeros(iter,1);
    selection=PSOMyself(tr_x, tr_y);
    str1={'FE', 'FS', 'NONE'};
    str2={'RBF1','SVM', 'RBF2'};%, 'KNN', 'DTree'
    for i=1:length(str1)
        for j=1:length(str2)
            str{(i-1)*length(str2)+j, 1}=str1{i};
            str{(i-1)*length(str2)+j, 2}=str2{j};
        end
    end

    %% 模型更新，不断找使EI最大的点
    for i=1:iter
        %% 选择训练数据
        N=size(tr_x,1);
        if N <=L
            disp(sprintf('No training data decrease'));
            tr_xx=tr_x;tr_yy=tr_y;
        else
            paixu = non_domination_sort_mod([tr_x, tr_y]);        
            data=paixu(1:floor(L/2), 1:V+M);
            paixu=paixu(floor(L/2)+1:end, 1:V+M);
            index=randperm(size(paixu,1));
            data=[data;paixu(index(1:L-floor(L/2)), 1:V+M)];
            tr_xx=data(:, 1:V);tr_yy=data(:, V+1:V+M);
        end
        %% 训练模型
        t1=clock;
        models=trainmodel(tr_xx, tr_yy, selection, str);
        timestore(i)=etime(clock, t1);
        %% 优化
        xnew=nsga_2(tr_xx, tr_yy, models, str);
        ynew=obj_real(xnew,Problem);
        %Ke=size(ynew,1);
        %ynew=ynew+repmat((max(tr_y)-min(tr_y)), Ke, 1).*(2*rand(Ke,M)-1)*0.1;
        tr_x=[tr_x;xnew];tr_y=[tr_y;ynew];
        disp(sprintf('%s :  %u / %u inner loop of %u / %u outside loop finished \n', Problem, i, iter, ceng, itr));
    end
    chrom(:,1:V)=tr_x;chrom(:,V+1:V+M)=tr_y;
