function [Ypred, PosCov, LCL, UCL]=BMGPPredictX_general(x,No,alphaopt,ropt,a0opt,s2,QI,QIRes,SD,n,k) 
%checked 3
m=size(x,1);
a=[x(:,2)*73+27 ones(m,1)*27];
Sx=x;%(x-repmat(Min,m,1))./repmat(Range,m,1)
Dist2=zeros(m,k);
for i=1:k        
    Dist2(:,i)=(0.5./(x(:,i)+0.5)-1).^2;   
end
w1=sum(Dist2,2);
IDist2=1./Dist2;
w2=sum(IDist2,2);
Index1=logical([zeros(m,1) Dist2==0]);

%unscaledsigmax=(prod((Dist2).^etaopt,2));
Q2=zeros(m,n);
for i=1:m
    for j=1:n
        Q2(i,j)=CompCorr(Sx(i,:),SD(j,:),ropt,No);
    end
end
w3=alphaopt*IDist2;
w=w1+alphaopt*w2;
Lambda=zeros(m,k+1);
Lambda(:,1)=w1./w;
Lambda(:,2:end)=w3./repmat(w,1,k);
Lambda(Index1)=1;
M=sum(Index1,2); IM=M>1;
if(sum(IM)>0)
Lambda(IM,Index1(IM,:))=1/M(IM);
end

miux=sum(Lambda(:,2:end).*a,2)+Lambda(:,1)*a0opt;
Ypred=miux+Q2*QIRes;
Q3=zeros(m,m);
for i=1:m
    for j=i:m
        Q3(i,j)=CompCorr(Sx(i,:),Sx(j,:),ropt,No);
        Q3(j,i)=Q3(i,j);        
    end
end
PosCov=s2*(Q3-Q2*QI*Q2');
RMSE=sqrt(max(diag(PosCov),0));
LCL=Ypred-norminv(0.99)*RMSE;
UCL=Ypred+norminv(0.99)*RMSE;
 

function Corr=CompCorr(x1,x2,r,No)
if(No==0)
    phi=r(3:4);
    r=r(1:2);
    Corr=prod(((abs(phi.*x1).^(r))+(abs(phi.*x2).^(r))-(abs(phi.*x1-phi.*x2).^(r)))/2);
elseif(No==1)
    Corr=prod(sinh(r.*(min(x1,x2))).*exp(-r.*(max(x1,x2))));
    
elseif(No==2)
    Corr=prod(min(x1,x2));
end