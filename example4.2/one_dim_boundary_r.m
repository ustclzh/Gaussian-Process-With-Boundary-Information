clear
clc
warning off
x=[0.1,0.13,0.16,0.2,0.3,0.5]';
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
len = 10.0; % side length for the square plate
hmax = len/20; % mesh size parameter
pres = 2; % external pressure
wMax_true=@(x) -0.0138*pres*len^4./(E*(x).^3);
wMax = wMax_true(x);
wMax_b=wMax_true(1);
theta_gp=GP_fit(x,wMax,wMax_b);
theta_bmgp=BMGP_fit(x,wMax,wMax_b);
theta_gp0=GPM_fit(x,wMax,wMax_b,0);
theta_gp1=GPM_fit(x,wMax,wMax_b,1);
theta_gp2=GPM_fit(x,wMax,wMax_b,2);
%%
x_test=min(x):0.01:1;
n=length(x_test);
y_p=zeros(n,5);
sig=zeros(n,5);
y_t=zeros(1,n);
for i=1:n
res=GP_predict(x,wMax,x_test(i),theta_gp,theta_bmgp,theta_gp0,theta_gp1,theta_gp2,wMax_b);
y_p(i,:)=res(1,:);
sig(i,:)=res(2,:);
y_t(i)=wMax_true(x_test(i));
end
p1=plot(x_test,y_p(:,1),'b','LineWidth',2); hold on
p5=plot(x_test,y_p(:,1)+(norminv(0.99)*sig(:,1)).^(1/2),'b--','LineWidth',2);
plot(x_test,y_p(:,1)-(norminv(0.99)*sig(:,1)).^(1/2),'b--','LineWidth',2)
p2=plot(x_test,y_p(:,2),'r','LineWidth',2);
p6=plot(x_test,y_p(:,2)+(norminv(0.99)*sig(:,2)).^(1/2),'r--','LineWidth',2);
plot(x_test,y_p(:,2)-(norminv(0.99)*sig(:,2)).^(1/2),'r--','LineWidth',2);
p3=plot(x,wMax,'ko','MarkerSize',12);
p4=plot(x_test,y_t,'kx','MarkerSize',12);
label={'Posterior mean of standard GP emulator','Posterior mean of BMGP emulator','Data','True function','LCL/UCL for standard GP emulator ','LCL/UCL for BMGP emulator '};
legend([p1,p2,p3,p4,p5,p6],label);hold on
ylabel('Deflection')
xlabel('Thickness')
print(gcf,'-dtiff','-r300','result_boundary_r');


function theta_gp=GP_fit(design,y,wMax_b)% all column
%%
global design_x y_out y_left y_right type_corr
type_corr=3;
y_left=y(1);
y_right=y(end);
design_x=[design;1];
y_out=[y;wMax_b];
a=(-5:1:5)';
problem = createOptimProblem('fminunc','x0',zeros(1,1),'objective',@lik_gp);
tpoints=CustomStartPointSet(a);
ms=MultiStart('StartPointsToRun','all','Display','off');
[theta_gp,fval0,exitflag0,solution0]=run(ms,problem,tpoints);
end

 
function theta_gp=GPM_fit(design,y,wMax_b,type)% all column
%%
global design_x y_out type_corr bound_r
type_corr=type;
design_x=design;
y_out=y;
bound_r=wMax_b;
switch type
    case 0
        a=sobolset(2);
        b=(net(a,10)-0.5);
        problem = createOptimProblem('fminunc','x0',zeros(1,2),'objective',@lik_gpmr);
        tpoints=CustomStartPointSet(b);
        ms=MultiStart('StartPointsToRun','all','Display','off');
        [theta_gp,fval0,exitflag0,solution0]=run(ms,problem,tpoints);
    case 1
        a=(-5:1:5)';
        problem = createOptimProblem('fminunc','x0',zeros(1,1),'objective',@lik_gpmr);
        tpoints=CustomStartPointSet(a);
        ms=MultiStart('StartPointsToRun','all','Display','off');
        [theta_gp,fval0,exitflag0,solution0]=run(ms,problem,tpoints);

    case 2
        theta_gp=1;
end
end




function y=GP_predict(design,y,x_p,theta_gp,theta_bmgp,theta_gp0,theta_gp1,theta_gp2,wMax_b)
theta_bmgp(2)=0;
[y1,sig1]=predict_gp(design,x_p,y,theta_gp,wMax_b);%column
[y2,sig2]=predict_bmgp(design,x_p,y,theta_bmgp,wMax_b);%column
[y3,sig3]=predict_gpm(design,x_p,y,theta_gp0,wMax_b,0);%column
[y4,sig4]=predict_gpm(design,x_p,y,theta_gp1,wMax_b,1);%column
[y5,sig5]=predict_gpm(design,x_p,y,theta_gp2,wMax_b,2);%column
y=[y1,sig1;y2,sig2;y3,sig3;y4,sig4;y5,sig5]';
end

function y=lik_gp(theta)
global design_x y_out
if (exp(theta)<0.1)||(exp(theta)>50)
    y=inf;
    return 
end
n=size(design_x,1);
R=corr_m(theta);
mu=sum(R\(y_out))/(sum(R\ones(n,1)));
sig=(y_out-mu)'*(R\(y_out-mu))/n;
y=n*log(sig)+log(det(R));
end


function y=lik_gpmr(theta)
% likelihood with right bound
global design_x y_out bound_r type_corr
design_x=1-design_x;
if type_corr==0
    if (exp(theta(1))>=2)
        y=inf;
        return
    end
else
    if (exp(theta(1))<0.1)||(exp(theta(1))>50)
        y=inf;
        return 
    end
end
n=size(design_x,1);
R=corr_m(theta);
sig=(y_out-bound_r)'*(R\(y_out-bound_r))/n;
y=n*log(sig)+log(det(R));
end

function theta_bmgp=BMGP_fit(design,y,wMax_b)% all column
%%
global design_x y_out y_left y_right wMax_r
wMax_r=wMax_b;
y_left=y(1);
y_right=y(end);
design_x=design;
y_out=y;
b=sobolset(3);
b=(net(b,30)-0.5);
problem = createOptimProblem('fminunc','x0',zeros(1,3),'objective',@lik_bmgp);
tpoints=CustomStartPointSet(b);
ms=MultiStart('StartPointsToRun','all','Display','off');
[theta_bmgp,fval1,exitflag1,solution1]=run(ms,problem,tpoints);
end
function y=lik_bmgp(theta)
global design_x y_out wMax_r
d=size(design_x,2);
n=size(y_out,2);
if (exp(theta(1))<0.1)||(exp(theta(1))>50)
    y=inf;
    return 
end
Q=corr_inf(theta);
D=abs(1./abs(1+2*design_x)-1/3);
alpha =theta(d+2);
lambda=D.^2./(D.^2+exp(alpha)./D.^2);
lambda1=1-lambda;
mu=(lambda'*(Q\lambda))\(lambda'*(Q\(y_out-lambda1*wMax_r)));
sig=(y_out-mu*lambda-lambda1*wMax_r)'*(Q\(y_out-mu*lambda-lambda1*wMax_r));
y=n*log(sig)+log(det(Q));
end
function y=corr_inf(theta)%specified
global design_x
d=size(design_x,2);
eta=theta(d+1);
D=abs(1./abs(1+2*design_x)-1/3);
D=D*D';
y=D.^exp(eta).*corr_m(theta(1:d));
end
function y=corr_m(theta)
global design_x type_corr
design=design_x;
n=size(design,1);
R=ones(n,n);
for i=1:n
    for j=1:n
        R(i,j)=corr_custom(design(i,:),design(j,:),theta,type_corr);
        R(j,i)=R(i,j);
    end
end
y=R;
end
function [y,sig]=predict_gp(design,x_p,y,theta,wMax_b)%column
global design_x y_out type_corr
type_corr=3;
design_x=[design;1];
n=size(design_x,1);
y_out=[y;wMax_b];
R=corr_m(theta);
mu=sum(R\(y_out))/(sum(sum(inv(R))));
r=ones(n,1);
for i=1:n
   r(i)=corr_custom(design_x(i,:),x_p,theta,3);
end
sig=max(0,(1-r'*(R\r))*(y_out-mu)'*(R\(y_out-mu))/n);
y=mu+r'*(R\(y_out-mu));
end
function [y,sig]=predict_gpm(design,x_p,y,theta,wMax_b,type)%column
global design_x y_out type_corr
type_corr=type;
design_x=1-design;
x_p=1-x_p;
n=size(design_x,1);
y_out=y;
R=corr_m(theta)+0.00000001*eye(n);
mu=wMax_b;
r=ones(n,1);
for i=1:n
   r(i)=corr_custom(design_x(i,:),x_p,theta,type_corr);
end
sig=((corr_custom(x_p,x_p,theta,type_corr)-r'*(R\r)))*(y_out-mu)'*(R\(y_out-mu))/n;
y=mu+r'*(R\(y_out-mu));
end
function [y,sig]=predict_bmgp(design,x_p,y,theta,wMax_b)%column
global design_x y_out type_corr
type_corr=3;
design_x=design;
y_out=y;
d=size(design,2);
eta=theta(d+1);
n=size(design,1);
Q=corr_inf(theta);
r=ones(n,1);
for i=1:n
   r(i)=corr_custom(design(i,:),x_p,theta(1:d),3);
end

D=abs(1./abs(1+2*design)-1/3);

alpha =theta(d+2);
lambda=D.^2./(D.^2+exp(alpha)./D.^2);
lambda1=1-lambda;
mu=(lambda'*(Q\lambda))\(lambda'*(Q\(y_out-lambda1*wMax_b)));
d=abs(1./abs(1+2*x_p)-1/3);
q=r.*(d*D).^exp(eta);
sig=(d^(2*exp(eta))-q'*(Q\q))*(y_out-mu*lambda-lambda1*wMax_b)'*(Q\(y_out-mu*lambda-lambda1*wMax_b))/n;
y=mu*d^2/(d^2+exp(alpha)/d^2)+(1-d^2/(d^2+exp(alpha)/d^2))*wMax_b+q'*(Q\(y_out-mu*lambda-wMax_b*lambda1));
end

function y=corr_custom(x1,x2,r,type)
%type: type of correlation function 
%
%
%
switch type 
    case 0
        d=size(x1,2);
        phi=exp(r(d+1:end));
        r=exp(r(1:d));
        y=prod(((abs(phi.*x1).^(r))+(abs(phi.*x2).^(r))-(abs(phi.*x1-phi.*x2).^(r)))/2);
    case 1
        y=prod(sinh(exp(r).*(min(x1,x2))).*exp(-exp(r).*(max(x1,x2))));
    case 2
        y=prod(min(x1,x2));
    case 3
        y=exp(-abs(x1-x2)./exp(r)).*(1+abs(x1-x2)./exp(r));
end

end
