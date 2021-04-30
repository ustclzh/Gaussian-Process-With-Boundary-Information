function [ropt0,RI0,beta0,RIRes0,sigma20,SD0out,n0out,Minout,Rangeout]=GPFit_general(Din,Yin,Noin)
%checked 3
global SD0 Y0 F n0 No0 Min Range gamma Temp_range
gamma=2;
No0=Noin;
D=Din;
if Noin==1
mu=rbf(D)';
else
mu=Din(:,2)*73+27;
end
Y0=Yin-mu;
[n0, d]=size(D); n0out=n0;
Min=min(D,[],1); Minout=Min;
Range=range(D,1); Rangeout=Range;
SD0=D;
SD0out=SD0;
F=ones(n0,1);

options=optimoptions(@fminunc,'MaxIter',10^5,'TolX',10^-6,'TolFun',10^-8,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton');

if(No0==0)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(0.91*ones(1,d+2))],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(1.9*ones(1,d+2))],options);      
    if((exitflag1<=0)||(exitflag2<=0))
        display('Likelihood optimization failure: algorithm did not converge.')   
    end
    ropt0=(fval1<fval2)*exp(paropt1)+(fval1>=fval2)*exp(paropt2);
elseif(No0==1)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(0.91*ones(1,d))],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(2.55*ones(1,d))],options);     
    if((exitflag1<=0)||(exitflag2<=0))
        display('Likelihood optimization failure: algorithm did not converge.')   
    end
    ropt0=(fval1<fval2)*exp(paropt1)+(fval1>=fval2)*exp(paropt2);
end
if (No0==2)
    ropt0=[];
else
disp('Optimal values of GP parameters')
disp(ropt0)
end
R=zeros(n0,n0);
for i=1:n0
    for j=i:n0
        R(i,j)=CompCorr(SD0(i,:),SD0(j,:),ropt0);
        R(j,i)=R(i,j);        
    end
end
RI0=invandlogdet(R);
RIF=RI0*F;
beta0=(RIF'*Y0)/(F'*RIF);
Res=Y0;%-F*beta0;
RIRes0=RI0*Res;
sigma20=Res'*RIRes0/n0;

function Objective=Obj(par)
global SD0 Y0 F n0 No0

if(sum(isnan(par))>0)
    Objective=Inf;
    return
end 
if(No0==0)
    r=exp(par);
    if((sum(r<0.01)>0)||(max(r(1:2))>=2))
        Objective=Inf;
        return
    end
elseif(No0==1)
    r=exp(par); 
    if((sum(r<0.01)>0)||(sum(r>50)>0))
        Objective=Inf;
        return        
    end      
end
R=zeros(n0,n0);
for i=1:n0
    for j=i:n0
        R(i,j)=CompCorr(SD0(i,:),SD0(j,:),r);
        R(j,i)=R(i,j);
    end
end
[RI0, LDR]=invandlogdet(R);
%RIF=RI0*F;
%beta0=(RIF'*Y0)/(F'*RIF);
Res=Y0;%-F*beta0;
sigma20=Res'*(RI0*Res)/(n0);
Objective=n0*log(max(sigma20,0))+LDR;

function Corr=CompCorr(x1,x2,r)
global No0 
if No0==0
    phi=r(3:4);
    r=r(1:2);
    Corr=prod(((abs(phi.*x1).^(r))+(abs(phi.*x2).^(r))-(abs(phi.*x1-phi.*x2).^(r)))/2);
end
if No0==1
    Corr=prod(sinh(r.*(min(x1,x2))).*exp(-r.*(max(x1,x2))));
end
if No0==2
    Corr=prod(min(x1,x2));
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
function y=rb(x1,x2)
global gamma
y=max((1-(sum((x1-x2).^2))^(1/2)),0)^gamma;





