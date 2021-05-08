clear all
%%
addpath(genpath(pwd));
load('design.mat')
for exp_num=1:60
exp_num
global hmax Temp_range halfl radius Ytest alpha_cl
hmax=30;
alpha_cl=0.99;
Temp_range=700;
halfl=0.45;
radius=0.075;
D=D_training{exp_num};
n=size(D,1);
Duncoded=[D(:,1)*hmax D(:,2)*Temp_range+20];
Y=zeros(n,1);
for i=1:n
    Y(i)=rodexample(Duncoded(i,1),Duncoded(i,2));
end
Dtest=D_Testing{exp_num};
%Dtest=[lhsdesign(40,2)];
ntest=size(Dtest,1);
Dtestuncoded=[Dtest(:,1)*hmax Dtest(:,2)*Temp_range+20];
Ytest=zeros(ntest,1);
for i=1:ntest
    Ytest(i)=rodexample(Dtestuncoded(i,1),Dtestuncoded(i,2));
end
%
[ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range,Condnum_std]=GPFit0(D,Y,1);
[Ypred0, PosCov0, LCL0, UCL0]=GPPredictX0(Dtest,1,ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range);

[alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k,Condnum_bm]=BMGPFit_ls(D,Y,1);
[Ypred, PosCov, LCL, UCL]=BMGPPredictX(Dtest,1,alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k);

[alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k]=FBGPFit(D,Y,0);
[YpredBmu, PosCovBmu, LCLBmu, UCLBmu]=FBGPPredictX(Dtest,0,alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k);

[roptDingmu,RIDingmu,betaDingmu,RIResDingmu,sigma2Dingmu,SDDingmu,nDingmu,MinDingmu,RangeDingmu]=BdryGPFit(D,Y,1);
[YpredDingmu, PosCovDingmu, LCLDingmu, UCLDingmu]=BdryGPPredictX(Dtest,1,roptDingmu,RIDingmu,betaDingmu,RIResDingmu,sigma2Dingmu,SDDingmu,nDingmu,MinDingmu,RangeDingmu);

[roptGraepel,RIGraepel,RIResGraepel,sigma2Graepel,SDGraepel,nGraepel]=GGPFit(D,Y,1);
[YpredGraepel, PosCovGraepel, LCLGraepel, UCLGraepel]=GGPPredictX(Dtest,1,roptGraepel,RIGraepel,RIResGraepel,sigma2Graepel,SDGraepel,nGraepel);
%
Condnum(exp_num,:)=[Condnum_std,Condnum_bm];
MAE0=mean(abs(Ytest-Ypred0));
RMSE0=sqrt(mean((Ytest-Ypred0).^2));
ALPI0=mean(UCL0-LCL0);
Coverage0=mean((Ytest>=LCL0).*(Ytest<=UCL0));

MAE=mean(abs(Ytest-Ypred));
RMSE=sqrt(mean((Ytest-Ypred).^2));
ALPI=mean(UCL-LCL);
Coverage=mean((Ytest>=LCL).*(Ytest<=UCL));
 
MAEBmu=mean(abs(Ytest-YpredBmu));
RMSEBmu=sqrt(mean((Ytest-YpredBmu).^2));
ALPIBmu=mean(UCLBmu-LCLBmu);
CoverageBmu=mean((Ytest>=LCLBmu).*(Ytest<=UCLBmu));

MAEDingmu=mean(abs(Ytest-YpredDingmu));
RMSEDingmuu=sqrt(mean((Ytest-YpredDingmu).^2));
ALPIDingmuu=mean(UCLDingmu-LCLDingmu);
CoverageDingmu=mean((Ytest>=LCLDingmu).*(Ytest<=UCLDingmu));

MAEGraepel=mean(abs(Ytest-YpredGraepel));
RMSEGraepel=sqrt(mean((Ytest-YpredGraepel).^2));
ALPIGraepel=mean(UCLGraepel-LCLGraepel);
CoverageGraepel=mean((Ytest>=LCLGraepel).*(Ytest<=UCLGraepel));
%
res=[MAE0  MAEBmu  MAEDingmu  MAE  MAEGraepel  
RMSE0  RMSEBmu  RMSEDingmuu RMSE  RMSEGraepel 
ALPI0  ALPIBmu  ALPIDingmuu  ALPI ALPIGraepel 
Coverage0  CoverageBmu  CoverageDingmu  Coverage CoverageGraepel ]';
RES{exp_num}=res;
Y_training(:,exp_num)=Y;
Y_testing(:,exp_num)=Ytest;
end
%%
for i=1:60
    A(i,:)=RES{i}(:)';
end
res.summary=reshape(mean(A),5,4);
%stdd=reshape(std(A),5,4)
%stdd=reshape(std(A(31:60,:))/sqrt(30),5,4)
res.se=reshape(std(A)/sqrt(60),5,4);
save('result.mat','res')




