function [Ypred, PosCov, LCL, UCL]=BMGPPredictX_Graepel(x,No,alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k) 
%checked 3
global Temp_range
m=size(x,1);
a=[x(:,2)*Temp_range+27 ones(m,1)*27];
Sx=(x-repmat(Min,m,1))./repmat(Range,m,1);
Dist2=zeros(m,k);
for i=1:k        
    Dist2(:,i)=(0.5./(x(:,i)+0.5)-1).^2;   
end
w1=sum(Dist2,2);
IDist2=1./Dist2;
w2=sum(IDist2,2);
Index1=logical([zeros(m,1) Dist2==0]);

unscaledsigmax=(prod((Dist2).^etaopt,2));
Q2=zeros(m,n);
for i=1:m
    for j=1:n
        Q2(i,j)=unscaledsigmax(i)*unscaledsigmaD(j)*CompCorr(Sx(i,:),SD(j,:),ropt,No);
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

mu=(1-x(:,1)).*a(:,1)+(1-x(:,2)).*(a(:,2)-(1-x(:,1))*27);

%miux=sum(Lambda(:,2:end).*a,2)+Lambda(:,1)*a0opt; 
Ypred=mu+Q2*QIRes;
if(m<=1000)
    Q3=zeros(m,m);
    for i=1:m
    for j=i:m
        Q3(i,j)=unscaledsigmax(i)*unscaledsigmax(j)*CompCorr(Sx(i,:),Sx(j,:),ropt,No);
        Q3(j,i)=Q3(i,j);        
    end
    end
    PosCov=s2*(Q3-Q2*QI*Q2');
    RMSE=sqrt(max(diag(PosCov),0));
else
    Q3=zeros(m,1);
    for i=1:m
        Q3(i)=(unscaledsigmax(i)^2)*CompCorr(Sx(i,:),Sx(i,:),ropt,No);
    end
    Q2T=Q2';
    Q4=zeros(m,1);
    for i=1:m
        Q4(i)=Q2(i,:)*(QI*Q2T(:,i));
    end      
   PosCov=s2*(Q3-Q4);
   RMSE=sqrt(max(PosCov,0));    
end
LCL=Ypred-norminv(0.99)*RMSE;
UCL=Ypred+norminv(0.99)*RMSE;

function Corr=CompCorr(x1,x2,r,No)
if(No==0)
    Corr=prod(r.^(abs(x1-x2)));
elseif(No==1)
    rho=x1-x2;
    rho1=sqrt(6)*abs(rho)./r;
    Corr=prod((exp(-rho1)).*(rho1+1));
elseif(No==2)
    rho=x1-x2;
    rho1=2*sqrt(2.5)*abs(rho)./r;
    Corr=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3),2);
end