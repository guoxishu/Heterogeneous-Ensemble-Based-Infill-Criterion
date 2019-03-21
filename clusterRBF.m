function [Centers Spreads W2 B2]=clusterRBF(SamIn, SamOut, ClusterNum)
    
    v = size(SamIn,2);
    Overlap = 1.0;% 隐节点重叠系数
    SamNum = size(SamIn,1);% 总样本数
    
    nn=0;
    while (nn==100||nn==0)
        index=randi([1,SamNum],ClusterNum,1);
        Centers = SamIn(index,:);%初始化聚类中心
        n=1;
        while n<100,
            NumberInClusters = zeros(ClusterNum,1); % 各类中的样本数，初始化为零
            IndexInClusters = zeros(ClusterNum,SamNum); % 各类所含样本的索引号
            % 按最小距离原则对所有样本进行分类
            for i = 1:SamNum
                AllDistance = dist(Centers,SamIn(i,:)');%计算第i个训练输入与所有聚类中心的距离
                [MinDist,Pos] = min(AllDistance);   %最小距离，该训练输入所属的聚类中心索引
                NumberInClusters(Pos) = NumberInClusters(Pos) + 1;
                IndexInClusters(Pos,NumberInClusters(Pos)) = i;%将属于该类的训练输入索引依次存入
            end   
            % 保存旧的聚类中心
            OldCenters = Centers;
            %重新计算聚类中心
            for i = 1:ClusterNum
                Index = IndexInClusters(i,1:NumberInClusters(i));%提取属于该类的训练输入索引
                Centers(i,:) = mean(SamIn(Index,:),1);    %取各类的平均作为新的聚类中心
            end
            % 判断新旧聚类中心是否一致，是则结束聚类
            EqualNum = sum(sum(Centers==OldCenters));%Centers与OldCenters所有对应位相减之后求总和
            if EqualNum ==v*ClusterNum,%表明新旧聚类中心一致
                break,
            end
            n=n+1;
        end
        nn=n;
    end
    % 计算各隐节点的扩展常数（宽度）
    AllDistances = dist(Centers,Centers'); % 计算隐节点数据中心间的距离（ClusterNum维数的方阵，对称矩阵）
    Maximum = max(max(AllDistances)); % 找出其中最大的一个距离
    for i = 1:ClusterNum % 将对角线上的0 替换为较大的值
        AllDistances(i,i) = Maximum+1;
    end
    Spreads = Overlap*min(AllDistances)'; % 以隐节点间的最小距离作为扩展常数; 并转换为列向量

    % 计算各隐节点的输出权值
    Distance = dist(Centers,SamIn'); % 计算各样本输入离各数据中心的距离(是ClusterNum X SamNum矩阵)
    SpreadsMat = repmat(Spreads,1,SamNum);%是ClusterNum X SamNum矩阵
    HiddenUnitOut = radbas(Distance./SpreadsMat); %计算隐节点输出阵; radbas是径向基传递函数
    HiddenUnitOutEx = [HiddenUnitOut' ones(SamNum,1)]'; % 考虑偏移(阈值)
    W2Ex = SamOut'*pinv(HiddenUnitOutEx); % 求广义输出权值
    W2 = W2Ex(:,1:ClusterNum); % 输出权值
    B2 = W2Ex(:,ClusterNum+1); % 偏移（阈值）
end