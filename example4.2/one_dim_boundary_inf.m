clear;clc;warning off; format long
x=[0.1,0.13,0.16,0.2,0.3,0.5]';
x=[0.4,0.5,0.6,0.7,0.8,0.9]';
x=[0.1,0.125,0.15,0.175,0.2]';
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
len = 10.0; % side length for the square plate
hmax = len/20; % mesh size parameter
pres = 2; % external pressure
wMax = -0.0138*pres*len^4./(E*(x).^3);
[theta_gp,theta_bmgp]=GP_fit(x,wMax);
x_test=min(x):0.01:max(x)+0.2;
n=length(x_test);
for i=1:n
res=GP_predict(x,wMax,x_test(i),theta_gp,theta_bmgp);
y_p(i,:)=res(1,:);
sig(i,:)=res(2,:);
y_t(i)=-0.0138*pres*len^4./(E*(x_test(i)).^3);
end
%%
p2=plot(x_test,y_p(:,2),'r','LineWidth',2); hold on
p6=plot(x_test,y_p(:,2)+(norminv(0.99)*sig(:,2)).^(1/2),'r','LineWidth',2);
plot(x_test,y_p(:,2)-(norminv(0.99)*sig(:,2)).^(1/2),'r','LineWidth',2);
p1=plot(x_test,y_p(:,1),'b--','LineWidth',2); hold on
p5=plot(x_test,y_p(:,1)+(norminv(0.99)*sig(:,1)).^(1/2),'b--','LineWidth',2);
plot(x_test,y_p(:,1)-(norminv(0.99)*sig(:,1)).^(1/2),'b--','LineWidth',2)
p3=plot(x,wMax,'ko','MarkerSize',12);
p4=plot(x_test,y_t,'kx','MarkerSize',12);
label={'Posterior mean / LCL / UCL of standard GP emulator','Posterior mean / LCL / UCL of BMGP emulator','Data','True function'};
legend([p1,p2,p3,p4],label,'FontSize',13,'Fontname', 'Times New Roman');hold on
ylabel('Maximum Deflection ','FontSize',13,'Fontname', 'Times New Roman')
xlabel('Thickness','FontSize',13,'Fontname', 'Times New Roman')
%xticklabels('FontSize',13)
a=2.2;
set(gca,'FontSize',13);
set(gcf,'Position',[200,200,1000,400]);
print(gcf,'-dtiff','-r300','result_boundary_inf');
%ylim([-0.5,0.28])
function [theta_gp,theta_bmgp]=GP_fit(design,y)% all column
%%
global design_x y_out alpha_opt
design_x=design;
y_out=y;
a=(-5:1:5)';
b=sobolset(2);
b=(net(b,20)-0.5);
problem = createOptimProblem('fminunc','x0',zeros(1,1),'objective',@lik_gp);
tpoints=CustomStartPointSet(a);
ms=MultiStart('StartPointsToRun','all','Display','off');
[theta_gp,fval0,exitflag0,solution0]=run(ms,problem,tpoints);


problem = createOptimProblem('fminunc','x0',zeros(1,1),'objective',@ls_bmgp_inf);
tpoints=CustomStartPointSet(a);
ms=MultiStart('StartPointsToRun','all','Display','off');
[alpha_opt,fval0,exitflag0,solution0]=run(ms,problem,tpoints);
problem = createOptimProblem('fminunc','x0',zeros(1,2),'objective',@lik_bmgp_inf);
tpoints=CustomStartPointSet(b);
ms=MultiStart('StartPointsToRun','all','Display','off');
[theta_bmgp,fval1,exitflag1,solution1]=run(ms,problem,tpoints);
end

function y=GP_predict(design,y,x_p,theta_gp,theta_bmgp)
global design_x y_out
design_x=design;
y_out=y;
[y1,sig1]=predict_gp(design,x_p,y,theta_gp);%column
[y2,sig2]=predict_bmgp(design,x_p,y,theta_bmgp);%column
y=[y1,sig1;y2,sig2]';
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
[invR,logdetR]=invandlogdet(R);
sig=(y_out-mu)'*(invR*(y_out-mu))/n;
y=n*log(sig)+logdetR;
end

function y=lik_bmgp_inf(theta)
global design_x y_out alpha_opt
n=size(design_x,1);
d=size(design_x,2);
if (exp(theta(1))<0.1)||(exp(theta(1))>50)
    y=inf;
    return 
end
Q=corr_inf(theta);
D=1./abs(1+2*design_x);
alpha =alpha_opt;
lambda=D.^2./(D.^2+exp(alpha)./D.^2);
%lamnda_1=1-lambda;
mu=(lambda'*(eye(n)\lambda))\(lambda'*(eye(n)\(y_out)));
[invQ,logdetQ]=invandlogdet(Q);
sig=(y_out-mu*lambda)'*(invQ*(y_out-mu*lambda))/n;
y=n*log(sig)+logdetQ;
end

function y=ls_bmgp_inf(alpha)
global design_x y_out
n=size(design_x,1);
D=1./abs(1+2*design_x);
lambda=D.^2./(D.^2+exp(alpha)./D.^2);
mu=(lambda'*(eye(n)\lambda))\(lambda'*(eye(n)\(y_out)));

y=sum((y_out-mu*lambda).^2);
end




function y=corr_m(theta)
global design_x
design=design_x;
d=size(design,2);
n=size(design,1);
R=ones(n,n);
for i=1:d
   R=R.*exp(-abs(design(:,i)- design(:,i)')/exp(theta(i))).*(1+abs(design(:,i)- design(:,i)')/exp(theta(i)));
end
y=R;
end

function y=corr_inf(theta)%specified
global design_x
d=size(design_x,2);
eta=theta(d+1);
D=1./abs(1+2*design_x);
D=D*D';
y=D.^exp(eta).*corr_m(theta(1:d));
end

function [y,sig]=predict_gp(design,x_p,y_out,theta)%column
d=size(design,2);
n=size(design,1);
R=corr_m(theta);
mu=sum(R\(y_out))/(sum(sum(inv(R))));
r=ones(n,1);
for i=1:d
   r=r.*exp(-abs(design(:,i)- x_p(i))/exp(theta(i))).*(1+abs(design(:,i)- x_p(i))/exp(theta(i)));
end
sig=max(0,(1-r'*(R\r))*(y_out-mu)'*(R\(y_out-mu))/n);
y=mu+r'*(R\(y_out-mu));
end

function [y,sig]=predict_bmgp(design,x_p,y_out,theta)%column
global alpha_opt

d=size(design,2);
eta=theta(d+1);
n=size(design,1);
Q=corr_inf(theta);
r=ones(n,1);
for i=1:d
   r=r.*exp(-abs(design(:,i)- x_p(i))/exp(theta(i))).*(1+abs(design(:,i)- x_p(i))/exp(theta(i)));
end
D=1./abs(1+2*design);
alpha =alpha_opt;
lambda=D.^2./(D.^2+exp(alpha)./D.^2);
mu=(lambda'*(eye(n)\lambda))\(lambda'*(eye(n)\(y_out)));
d=1./abs(1+2*x_p);
q=r.*(d*D).^exp(eta);
sig=max(0,(d^(2*exp(eta))-q'*(Q\q))*(y_out-mu*lambda)'*(Q\(y_out-mu*lambda))/n);
y=mu*d^2/(d^2+exp(alpha)/d^2)+q'*(Q\(y_out-mu*lambda));
end

function [invR, logdetR]=invandlogdet(R)
%checke
%Inputs: R=positive definite matrix. 
%Outputs: invR=inverse of R, logdetR=natural logarithm of the determinant of R.
[CR, p]=chol(R);
if(p~=0)
    [Eigvec, Eigval]=svd(R);
    Eigval=diag(Eigval);
    invR=Eigvec*diag(1./Eigval)*Eigvec';
    logdetR=sum(log(Eigval));
else
    ICR=CR\speye(size(CR,1));
    invR=ICR*(ICR');
    logdetR=2*sum(log(diag(CR)));
end
end
