function [Ypred0, PosCov0, LCL0, UCL0]=BdryGPPredictX(x,No,ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range)
m=size(x,1);
global Ytest Temp_range
ropt0=ropt0(1:2);
Sx=x;
R2=zeros(m,n0);
for i=1:m
    for j=1:n0
        R2(i,j)=CompCorr(Sx(i,:),SD0(j,:),ropt0,No);
    end
end
if No==1
    mu=rbf(x)';
else
    mu=x(:,2)*73+27;
end
%
%mu-Ytest
Ypred0=R2*RIRes0+mu;
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
end

function Corr=CompCorr(x1,x2,r,No)
if No==0
    phi=r(3:4);
    r=r(1:2);
    Corr=prod(((abs(phi.*x1).^(r))+(abs(phi.*x2).^(r))-(abs(phi.*x1-phi.*x2).^(r)))/2);
end
if No==1
    Corr=prod(sinh(r.*(min(x1,x2))).*exp(-r.*(max(x1,x2))));
end
if No==2
    Corr=prod(min(x1,x2));
end

end
function y=rbf(D)
global Temp_range
D1=D(:,1);
D2=D(:,2);
n=size(D,1);
for i=1:n
    p=[D2(i)*Temp_range+20,20];
    phi=[rb(D1(i),0),rb(D2(i),0)];
    Phi=[rb(0,0),rb([D1(i),0],[0,D2(i)]);rb([D1(i),0],[0,D2(i)]),rb(0,0)];
    y(i)=phi*(Phi\p');
end
end
function y=rb(x1,x2)
global gamma
y=max((1-(sum((x1-x2).^2))^(1/2)),0)^(gamma);
end