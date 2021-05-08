function [ropt0,RI0,beta0,RIRes0,sigma20,SD0out,n0out,Minout,Rangeout,Condnum]=GPFit0(Din,Yin,Noin)
%checked 3
global SD0 Y0 F n0 No0 Min Range

No0=Noin;
D=Din;
Y0=Yin;
[n0, d]=size(D); n0out=n0;
Min=min(D,[],1); Minout=Min;
Range=range(D,1); Rangeout=Range;
SD0=(D-repmat(Min,n0,1))./repmat(Range,n0,1); SD0out=SD0;
F=ones(n0,1);

options=optimoptions(@fminunc,'MaxIter',10^5,'TolX',10^-6,'TolFun',10^-8,'MaxFunEvals',10^5,'Display','off','Algorithm','quasi-newton');

if(No0==0)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[norminv(0.25*ones(1,d))],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[norminv(0.75*ones(1,d))],options);
elseif(No0==1)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,[log(0.91*ones(1,d))],options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,[log(2.55*ones(1,d))],options);      
elseif(No0==2)
    [paropt1, fval1, exitflag1]=fminunc(@Obj,log(0.89*ones(1,d)),options); 
    [paropt2, fval2, exitflag2]=fminunc(@Obj,log(2.24*ones(1,d)),options);    
end

if((exitflag1<=0)||(exitflag2<=0))
    display('Likelihood optimization failure: algorithm did not converge.')   
end
if(fval1<fval2)
    if(No0==0)
        ropt0=normcdf(paropt1);
    elseif(No0==1)
        ropt0=exp(paropt1);        
    elseif(No0==2)
        ropt0=exp(paropt1);        
    end
else
    if(No0==0)
        ropt0=normcdf(paropt2);
    elseif(No0==1)
        ropt0=exp(paropt2);  
    elseif(No0==2)
        ropt0=exp(paropt2);          
    end
end
display('Optimal values of GP parameters')
disp(ropt0)
R=zeros(n0,n0);
for i=1:n0
    for j=i:n0
        R(i,j)=CompCorr(SD0(i,:),SD0(j,:),ropt0);
        R(j,i)=R(i,j);        
    end
end
RI0=invandlogdet(R);
Condnum=cond(R);
RIF=RI0*F;
beta0=(RIF'*Y0)/(F'*RIF);
Res=Y0-F*beta0;
RIRes0=RI0*Res;
sigma20=Res'*RIRes0/n0;

function Objective=Obj(par)
global SD0 Y0 F n0 No0

if(sum(isnan(par))>0)
    Objective=Inf;
    return
end 
if(No0==0)
    r=normcdf(par);
    if((sum(r<0.01)>0)||(sum(r>0.999)>0))
        Objective=Inf;
        return
    end
elseif(No0==1)
    r=exp(par);    
    if((sum(r<0.01)>0)||(sum(r>53.9511207457687)>0))
        Objective=Inf;
        return        
    end      
elseif(No0==2)
    r=exp(par);
    if((sum(r<0.291731110468193)>0)||(sum(r>40.7953930912641)>0))
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
RIF=RI0*F;
beta0=(RIF'*Y0)/(F'*RIF);
Res=Y0-F*beta0;
sigma20=Res'*(RI0*Res)/(n0);
Objective=n0*log(max(sigma20,0))+LDR;

function Corr=CompCorr(x1,x2,r)
global No0
if(No0==0)
    Corr=prod(r.^(abs(x1-x2)));
elseif(No0==1)
    rho=x1-x2;
    rho1=sqrt(6)*abs(rho)./r;
    Corr=prod((exp(-rho1)).*(rho1+1));
elseif(No0==2)
    rho=x1-x2;
    rho1=2*sqrt(2.5)*abs(rho)./r;
    Corr=prod((exp(-rho1)/3).*(rho1.^2+3*rho1+3),2); 
end