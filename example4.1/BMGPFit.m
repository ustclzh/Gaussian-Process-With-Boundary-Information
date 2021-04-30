function [alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Minout,Rangeout,SDout,nout,kout]=BMGPFit(Din,Yin,Noin)
%checked 3
global D Min Range SD Y n No a k Dist2 IDist2 w1 w2 Index1 Temp_range 

No=Noin;
D=Din;
Y=Yin;
[n, d]=size(D); nout=n;
k=2; kout=k;
a=[D(:,2)*Temp_range+20 ones(n,1)*20];

Min=min(D,[],1); Minout=Min;
Range=range(D,1); Rangeout=Range;
SD=(D-repmat(Min,n,1))./repmat(Range,n,1); SDout=SD;

Dist2=zeros(n,k);
for i=1:k
    Dist2(:,i)=((0.5./(D(:,i)+0.5)-1).^2);
end

w1=sum(Dist2,2);
IDist2=1./Dist2;
w2=sum(IDist2,2);
Index1=logical([zeros(n,1) Dist2==0]);

options=optimoptions(@fminunc,'MaxIter',10^5,'TolX',10^-6,'TolFun',10^-8,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton');
Sum1=min(w1);
Sum2=min(w2);
startalpha0=Sum1/Sum2;

bestfval=Inf;
startalphaetacandidates=[startalpha0 0.1; startalpha0 0.7; startalpha0 1; startalpha0 3];
for j=1:4
startalpha=startalphaetacandidates(j,1);
starteta=startalphaetacandidates(j,2);
logstarteta=log(starteta);
unscaledsigmaD=prod(Dist2.^starteta,2);
if(sum(unscaledsigmaD<=0)==0)
if(No==0)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(startalpha) logstarteta norminv(0.25)*ones(1,d)],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(startalpha) logstarteta norminv(0.75)*ones(1,d)],options);
elseif(No==1)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(startalpha) logstarteta log(0.91)*ones(1,d)],options);  
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(startalpha) logstarteta log(2.55)*ones(1,d)],options);      
elseif(No==2)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(startalpha) logstarteta log(0.89)*ones(1,d)],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(startalpha) logstarteta log(2.24)*ones(1,d)],options);    
end

if((exitflag1<=0)||(exitflag2<=0))
    display('Likelihood optimization failure: algorithm did not converge.')   
end
if(min(fval1,fval2)<bestfval)
if(fval1<fval2)
    bestfval=fval1;
    if(No==0)
        alphaopt=exp(paropt1(1));
        etaopt=exp(paropt1(2));
        ropt=normcdf(paropt1((3):end));
    elseif(No==1)
        alphaopt=exp(paropt1(1));
        etaopt=exp(paropt1(2));           
        ropt=exp(paropt1((3):end));        
    elseif(No==2)
        alphaopt=exp(paropt1(1));
        etaopt=exp(paropt1(2));          
        ropt=exp(paropt1((3):end));        
    end
else
    bestfval=fval2;
    if(No==0)
        alphaopt=exp(paropt2(1));
        etaopt=exp(paropt2(2));        
        ropt=normcdf(paropt2((3):end));
    elseif(No==1)
        alphaopt=exp(paropt2(1));
        etaopt=exp(paropt2(2));       
        ropt=exp(paropt2((3):end));  
    elseif(No==2)
        alphaopt=exp(paropt2(1));
        etaopt=exp(paropt2(2));       
        ropt=exp(paropt2((3):end));          
    end
end
end
end
end
display('Optimal values of correlation parameters')
disp(ropt)

unscaledsigmaD=(prod((Dist2).^etaopt,2));
R=zeros(n,n);
for i=1:n
    for j=i:n
        R(i,j)=CompCorr(SD(i,:),SD(j,:),ropt);
        R(j,i)=R(i,j);
    end
end
RI=invandlogdet(R);

sI=diag(1./unscaledsigmaD);
QI=sI*RI*sI;
%cond(invandlogdet(QI))
w3=alphaopt*IDist2;
w=w1+alphaopt*w2;
Lambda=zeros(n,k+1);
Lambda(:,1)=w1./w;
Lambda(:,2:end)=w3./repmat(w,1,k);
Lambda(Index1)=1;
M=sum(Index1,2); IM=M>1;
if(sum(IM)>0)
Lambda(IM,Index1(IM,:))=1/M(IM);
end
Z=Y-sum(Lambda(:,2:end).*a,2);
a0opt=(Lambda(:,1)'*QI*Z)/(Lambda(:,1)'*QI*Lambda(:,1));
Res=Z-Lambda(:,1)*a0opt;
QIRes=QI*Res;
s2=Res'*QIRes/n;

function Objective=Obj(par)
global SD Y n No a k Dist2 IDist2 w1 w2 Index1

if(sum(isnan(par))>0)
    Objective=Inf;
    return
end
alpha=exp(par(1));
eta=exp(par(2));
par=par((3):end);

if(No==0)
    r=normcdf(par);
    if((sum(r<0.001)>0)||(sum(r>0.999)>0))
        Objective=Inf;
        return
    end
elseif(No==1)
    r=exp(par);    
    if((sum(r<0.265285395160763)>0)||(sum(r>53.9511207457687)>0))
        Objective=Inf;
        return        
    end 
%     if eta<0.5
%         Objective=Inf;
%         return    
%     end
elseif(No==2)
    r=exp(par);
    if((sum(r<0.291731110468193)>0)||(sum(r>40.7953930912641)>0))
        Objective=Inf;
        return        
    end
end

unscaledsigmaD=(prod((Dist2).^eta,2));
R=zeros(n,n);
for i=1:n
    for j=i:n
        R(i,j)=CompCorr(SD(i,:),SD(j,:),r);
        R(j,i)=R(i,j);
    end
end

[RI, LDR]=invandlogdet(R);
sI=diag(1./unscaledsigmaD);
QI=sI*RI*sI;
LDR=2*sum(log(unscaledsigmaD))+LDR;

w3=alpha*IDist2;
w=w1+alpha*w2;
Lambda=zeros(n,k+1);
Lambda(:,1)=w1./w;
Lambda(:,2:end)=w3./repmat(w,1,k);
Lambda(Index1)=1;
M=sum(Index1,2); IM=M>1;
if(sum(IM)>0)
Lambda(IM,Index1(IM,:))=1/M(IM);
end
Z=Y-sum(Lambda(:,2:end).*a,2);
a0opt=(Lambda(:,1)'*QI*Z)/(Lambda(:,1)'*QI*Lambda(:,1));
Res=Z-Lambda(:,1)*a0opt;
s2=Res'*(QI*Res)/(n);
Objective=n*log(max(s2,0))+LDR;
    
function Corr=CompCorr(x1,x2,r)
global No
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