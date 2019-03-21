function y=SVMtest(x, p, N)

global M
ps1=p{1};ps2=p{2};
xtest0=mapminmax('apply',x',ps1);xtest0=xtest0';
for j=1:M
    for k=1:N
        py(k,j) = svmpredict(0,xtest0(k,:),p{j+2}, '-q');
    end
end
y=mapminmax('reverse',py',ps2);y=y';  