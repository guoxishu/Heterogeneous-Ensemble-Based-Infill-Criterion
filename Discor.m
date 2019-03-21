function D=Discor(x, y)

[n, x_n]=size(x);
%disp(sprintf('%u hang %u lie', n, x_n));
if ~isempty(y)
    a=dist(x,x');b=dist(y,y');
    A=a-repmat(mean(a,2),1,n)-repmat(mean(a,1),n,1)+mean(mean(a));
    B=b-repmat(mean(b,2),1,n)-repmat(mean(b,1),n,1)+mean(mean(b));
    covxy=sum(sum(A.*B))/(n^2);
    covx=sum(sum(A.^2))/(n^2);
    covy=sum(sum(B.^2))/(n^2);
    D=sqrt(covxy)/sqrt(sqrt(covx)*sqrt(covy));

else
    if x_n==1
        D=0;
    else
        for i=1:x_n
            xa=x(:,i);
            index=setdiff([1:x_n],i);
            xb=x(:,index);    
            a=dist(xa,xa');b=dist(xb,xb');
            A=a-repmat(mean(a,2),1,n)-repmat(mean(a,1),n,1)+mean(mean(a));
            B=b-repmat(mean(b,2),1,n)-repmat(mean(b,1),n,1)+mean(mean(b));
            covxy=sum(sum(A.*B))/(n^2);
            covx=sum(sum(A.^2))/(n^2);
            covy=sum(sum(B.^2))/(n^2);
            D(1,i)=sqrt(covxy)/sqrt(sqrt(covx)*sqrt(covy));
        end
        D=mean(D);
    end
         
end