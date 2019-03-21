function models=RBFMyself2(tr_x, tr_y)

global M V
%tic;
ClusterNum=round(sqrt(M+V)+3);
[models.Centers models.Spreads models.W2 models.B2] = ...
                                           clusterRBF(tr_x, tr_y,ClusterNum);
                                       %toc;
                                       %a