function f = obj_real(x, Problem)
    
    global M V
    k=6;l=V-k;%k必须是M-1的倍数；WFG2和WFG3，l必须是2的倍数
    switch Problem
        case {'DTLZ1'}   
            g = 100*(V+1-M+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            for i = 1 : M
                f(:,i) = 0.5.*prod(x(:,1:M-i),2).*(1+g);
                if i > 1
                    f(:,i) = f(:,i).*(1-x(:,M-i+1));
                end
            end       
        case 'DTLZ2'
            g = sum((x(:,M:end)-0.5).^2,2);
            for i = 1 : M
                f(:,i) = (1+g).*prod(cos(0.5.*pi.*x(:,1:M-i)),2);
                if i > 1
                    f(:,i) = f(:,i).*sin(0.5.*pi.*x(:,M-i+1));
                end
            end  
        case {'DTLZ3'}   
            g = 100*(V+1-M+sum((x(:,M:end)-0.5).^2-cos(20.*pi.*(x(:,M:end)-0.5)),2));
            for i = 1 : M
                f(:,i) = (1+g).*prod(cos(0.5.*pi.*x(:,1:M-i)),2);
                if i > 1
                    f(:,i) = f(:,i).*sin(0.5.*pi.*x(:,M-i+1));
                end
            end           
            
        case {'DTLZ4'} 
            alpha=100;
            g = sum((x(:,M:end)-0.5).^2,2);
            for i = 1 : M
                f(:,i) = (1+g).*prod(cos(0.5.*pi.*(x(:,1:M-i).^alpha)),2);
                if i > 1
                    f(:,i) = f(:,i).*sin(0.5.*pi.*(x(:,M-i+1).^alpha));
                end
            end              
            
        case {'DTLZ5'}  
            g = sum((x(:,M:end)-0.5).^2,2);
            x(:,2:M-1)=0.5*(1+2*repmat(g,1,M-2).*x(:,2:M-1))./(repmat(g,1,M-2)+1);
            for i = 1 : M
                f(:,i) = (1+g).*prod(cos(0.5.*pi.*x(:,1:M-i)),2);
                if i > 1
                    f(:,i) = f(:,i).*sin(0.5.*pi.*x(:,M-i+1));
                end
            end    
        
        case {'DTLZ6'}  
            g = sum(x(:,M:end).^0.1,2);
            x(:,2:M-1)=0.5*(1+2*repmat(g,1,M-2).*x(:,2:M-1))./(repmat(g,1,M-2)+1);
            for i = 1 : M
                f(:,i) = (1+g).*prod(cos(0.5.*pi.*x(:,1:M-i)),2);
                if i > 1
                    f(:,i) = f(:,i).*sin(0.5.*pi.*x(:,M-i+1));
                end
            end
            
        case {'DTLZ7'}  
            g = 1+9/(V+1-M)*sum(x(:,M:end),2);
            for i = 1 : M-1
                f(:,i) = x(:,i);
            end
            h=M-sum(f(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi*f(:,1:M-1))),2);
            f(:,M)=(1+g).*h;
            
        case {'ZDT6'}
            [m,n]=size(x);
            f(:,1) = 1-exp(-4*x(:,1)).*(sin(6*pi*x(:,1)).^6);
            summ = zeros(m,1);
            for i = 2 : n
                summ = summ + x(:,i);
            end
            g_x = 1 + 9*((summ)/(n-1)).^0.25;
            f(:,2) = g_x.*(1 - (f(:,1)./g_x).^2);   
            
        case {'ZDT3'}    
            f(:,1) = x(:,1);
            g_x = 1 +9*sum(x(:,2:end),2)/(size(x,2)-1);
            f(:,2) = g_x.*(1 - (x(:,1)./g_x).^0.5-x(:,1)./g_x.*sin(10*pi*x(:,1)));
            
        case {'ZDT4'}
            [m,n]=size(x);
            f(:,1) = x(:,1);
            summ = zeros(m,1);
            for i = 2 : n
                summ = summ + x(:,i).^2-10*cos(4*pi*x(:,i));
            end
            g_x = 1 +10*(n-1)+summ;
            f(:,2) = g_x.*(1 - (x(:,1)./g_x).^0.5);
        
        case {'WFG1'} 
            [N,V]=size(x);
            x=x./repmat(2*(1:V),N,1);
            %t1
            A=0.35;x(:,k+1:V)=abs(x(:,k+1:V)-A)./abs(floor(A-x(:,k+1:V))+A);
            %t2
            A=0.8;B=0.75;C=0.85;
            x(:,k+1:V)=A+min(0,floor(x(:,k+1:V)-B))*A/B.*(B-x(:,k+1:V))-...
                                   min(0,floor(C-x(:,k+1:V)))*(1-A)/(1-C).*(x(:,k+1:V)-C);
            %t3
            A=0.02;
            x=x.^A;
            %t4
            for i=1:M-1
                A=((i-1)*k/(M-1)+1:i*k/(M-1));
                t3(:,i)=sum(x(:, A).*repmat(2*A,size(x,1),1),2)/sum(2*A,2);
            end
            A=(k+1:V);
            t3(:,M)=sum(x(:,A).*repmat(2*A,size(x,1),1),2)/sum(2*A,2);
            
            alpha=1;A=5;
            for j=1:M-1
                h(:,j)=prod(1-cos(t3(:,1:M-j)*pi/2),2);   
                if j>1 
                    h(:,j)=prod(1-cos(t3(:,1:M-j)*pi/2),2).*(1-sin(t3(:,M-j+1)*pi/2));
                end
            end
            h(:,M)=(1-t3(:,1)-cos(2*A*pi*t3(:,1)+pi/2)/2/A/pi).^alpha;
            f=repmat(t3(:,M),1,M)+h.*repmat(2*(1:M),N,1);
            
        case {'WFG2'} 
            [N,V]=size(x);
            x=x./repmat(2*(1:V),N,1);
            %t1
            A=0.35;x(:,k+1:V)=abs(x(:,k+1:V)-A)./abs(floor(A-x(:,k+1:V))+A);
            %t2
            A=2;
            t2(:, 1:k)=x(:, 1:k);
            for i=1:l/2
                t2(:,k+i)=(x(:,k+2*i-1)+x(:,k+2*i)+2*abs(x(:,k+2*i-1)-x(:,k+2*i)))/3;
            end
            %t3
            for i=1:M-1
                t3(:,i)=mean(t2(:, (i-1)*k/(M-1)+1:i*k/(M-1)), 2);
            end
            t3(:,M)=mean(t2(:,k+1:k+l/2), 2);

            alpha=1;belta=1;A=5;
            for j=1:M
                h(:,j)=prod(1-cos(t3(:,1:M-j)*pi/2),2); 
                if j>1  
                    h(:,j)=prod(1-cos(t3(:,1:M-j)*pi/2),2).*(1-sin(t3(:,M-j+1)*pi/2));
                end
            end
            h(:,M)=1-t3(:,1).^alpha.*(cos(A*(t3(:,1)).^belta*pi)).^2;
            f=repmat(t3(:,M),1,M)+h.*repmat(2*(1:M),N,1);
            
        case {'WFG3'} 
            [N,V]=size(x);
            x=x./repmat(2*(1:V),N,1);
            %t1
            A=0.35;x(:,k+1:V)=abs(x(:,k+1:V)-A)./abs(floor(A-x(:,k+1:V))+A);
            %t2
            A=2;
            t2(:, 1:k)=x(:, 1:k);
            for i=1:l/2
                t2(:,k+i)=(x(:,k+2*i-1)+x(:,k+2*i)+2*abs(x(:,k+2*i-1)-x(:,k+2*i)))/3;
            end
            %t3
            for i=1:M-1
                t3(:,i)=mean(t2(:, (i-1)*k/(M-1)+1:i*k/(M-1)), 2);
            end
            t3(:,M)=mean(t2(:,k+1:k+l/2), 2);
            %x
            xx(:,1)=t3(:,1);
            for i=2:M-1
                xx(:,i)=t3(:,M).*(t3(:,i)-0.5)+0.5;
            end
            xx(:,M)=t3(:,M);
            
            for j=1:M
                h(:,j)=prod(xx(:,1:M-j),2);
                if j>1
                    h(:,j)=h(:,j).*(1-xx(:,M-j+1));
                end
            end
            f=repmat(xx(:,M),1,M)+h.*repmat(2*(1:M),N,1);   
            
        case {'WFG4'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);

            t1 = zeros(N,K+L);
            t1 = s_multi(z01,30,10,0.35);

            t2 = zeros(N,M);
            for i = 1 : M-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,M) = t2(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
        case {'WFG5'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1 = s_decept(z01,0.35,0.001,0.05);

            t2 = zeros(N,M);
            for i = 1 : M-1
                t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,M) = t2(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
        case {'WFG6'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1(:,1:K)     = z01(:,1:K);
            t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

            t2 = zeros(N,M);
            for i = 1 : M-1
                t2(:,i) = r_nonsep(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
            end
            % Same as <t2(:,M)=r_nonsep(t1(:,K+1:end),L)>
            SUM = zeros(N,1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t1(:,i)-t1(:,j));
                end
            end
            t2(:,M) = (sum(t1(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
            % -------------------------------------------

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
            end
            x(:,M) = t2(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
        case {'WFG7'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
            Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
            t1(:,1:K) = z01(:,1:K).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K)).*abs(floor(0.5-Y(:,1:K))+0.98/49.98)));
            % ------------------------------------------------------------------------------------------
            t1(:,K+1:end) = z01(:,K+1:end);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);

            t3 = zeros(N,M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,M) = t3(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
            
        case {'WFG8'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
            z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            t1(:,1:K) = z01(:,1:K);
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,1:i-1),ones(1,i-1)),0.98/49.98,0.02,50)>
            Y = (cumsum(z01,2)-z01)./repmat(0:K+L-1,N,1);
            t1(:,K+1:K+L) = z01(:,K+1:K+L).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,K+1:K+L)).*abs(floor(0.5-Y(:,K+1:K+L))+0.98/49.98))); 
            % --------------------------------------------------------------------------------------

            t2 = zeros(N,K+L);
            t2(:,1:K)     = t1(:,1:K);
            t2(:,K+1:end) = s_linear(t1(:,K+1:end),0.35);

            t3 = zeros(N,M);
            for i = 1 : M-1
                t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
            end
            t3(:,M) = r_sum(t2(:,K+1:K+L),ones(1,L));

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,M) = t3(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
     
        case {'WFG9'} 
            PopDec=x;
            N=size(PopDec,1);
            L=l;K=k;
            D = 1;
            S = 2 : 2 : 2*M;
            A = ones(1,M-1);
            
             z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
            
            t1 = zeros(N,K+L);
            % Same as <t1(:,i)=b_param(z01(:,i),r_sum(z01(:,i+1:end),ones(1,K+L-i)),0.98/49.98,0.02,50)>
            Y = (fliplr(cumsum(fliplr(z01),2))-z01)./repmat(K+L-1:-1:0,N,1);
            t1(:,1:K+L-1) = z01(:,1:K+L-1).^(0.02+(50-0.02)*(0.98/49.98-(1-2*Y(:,1:K+L-1)).*abs(floor(0.5-Y(:,1:K+L-1))+0.98/49.98)));
            % ------------------------------------------------------------------------------------------
            t1(:,end)     = z01(:,end);

            t2 = zeros(N,K+L);
            t2(:,1:K)     = s_decept(t1(:,1:K),0.35,0.001,0.05);
            t2(:,K+1:end) = s_multi(t1(:,K+1:end),30,95,0.35);

            t3 = zeros(N,M);
            for i = 1 : M-1
                t3(:,i) = r_nonsep(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),K/(M-1));
            end
            % Same as <t3(:,M)=r_nonsep(t2(:,K+1:end),L)>
            SUM = zeros(N,1);
            for i = K+1 : K+L-1
                for j = i+1 : K+L
                    SUM = SUM + abs(t2(:,i)-t2(:,j));
                end
            end
            t3(:,M) = (sum(t2(:,K+1:end),2)+SUM*2)/ceil(L/2)/(1+2*L-2*ceil(L/2));
            % -------------------------------------------

            x = zeros(N,M);
            for i = 1 : M-1
                x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
            end
            x(:,M) = t3(:,M);

            h = concave(x);
            f = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;
    end
end

function Output = s_multi(y,A,B,C)
    Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end

function Output = s_decept(y,A,B,C)
    Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end



