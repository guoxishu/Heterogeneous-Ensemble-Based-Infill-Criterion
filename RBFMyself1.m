function net=RBFMyself1(tr_x, tr_y)

%tic;
goal=sqrt(sum((max(tr_y)-min(tr_y)).^2))*0.05;
net=newrb(tr_x', tr_y', goal, 1.0, size(tr_x,2), 1e+06);%
%toc;
%a

%这两个是默认值，%spread,  Maximum number of neurons
%主要是为了设置后面的参数，不显示图% Number of neurons to add between displays