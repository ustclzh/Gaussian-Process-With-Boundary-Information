function [Ypred0, PosCov0, LCL0, UCL0]=GPPredictX0(x,No,ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range)
%checked 3

m=size(x,1);
Sx=(x-repmat(Min,m,1))./repmat(Range,m,1);
R2=zeros(m,n0);
for i=1:m
    for j=1:n0
        R2(i,j)=CompCorr(Sx(i,:),SD0(j,:),ropt0,No);
    end
end
Ypred0=beta0+R2*RIRes0;
if(m<=1000)
    R3=zeros(m,m);
    for i=1:m
    for j=i:m
        R3(i,j)=CompCorr(Sx(i,:),Sx(j,:),ropt0,No);
        R3(j,i)=R3(i,j);        
    end
    end
    PosCov0=sigma20*(R3-R2*RI0*R2');
    RMSE=sqrt(max(diag(PosCov0),0));
else
    R3=zeros(m,1);
    for i=1:m
        R3(i)=CompCorr(Sx(i,:),Sx(i,:),ropt0,No);
    end
    R2T=R2';
    R4=zeros(m,1);
    for i=1:m
        R4(i)=R2(i,:)*(RI0*R2T(:,i));
    end
    PosCov0=sigma20*(R3-R4);
    RMSE=sqrt(max(PosCov0,0));    
end
LCL0=Ypred0-norminv(0.99)*RMSE;
UCL0=Ypred0+norminv(0.99)*RMSE;

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