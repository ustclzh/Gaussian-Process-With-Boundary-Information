clear all
global hmax Temp_range halfl radius Ytest
hmax=30;
Temp_range=700;
halfl=0.45;
radius=0.1;
D=importdata('D1.csv');
n=size(D,1);
Duncoded=[D(:,1)*hmax D(:,2)*Temp_range+20];
Y=zeros(n,1);
for i=1:n
    Y(i)=rodexample(Duncoded(i,1),Duncoded(i,2));
end
Dtest=[importdata('D2.csv');importdata('D3.csv')];
ntest=size(Dtest,1);
Dtestuncoded=[Dtest(:,1)*hmax Dtest(:,2)*Temp_range+20];
Ytest=zeros(ntest,1);
for i=1:ntest
    Ytest(i)=rodexample(Dtestuncoded(i,1),Dtestuncoded(i,2));
end
[ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range]=GPFit0(D,Y,1);
[Ypred0, PosCov0, LCL0, UCL0]=GPPredictX0(Dtest,1,ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range);
%%
[alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k]=BMGPFit(D,Y,1);
[Ypred, PosCov, LCL, UCL]=BMGPPredictX(Dtest,1,alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k);

[alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k]=BMGPFit_general(D,Y,0);
[YpredBmu, PosCovBmu, LCLBmu, UCLBmu]=BMGPPredictX_general(Dtest,0,alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k);

[roptDingmu,RIDingmu,betaDingmu,RIResDingmu,sigma2Dingmu,SDDingmu,nDingmu,MinDingmu,RangeDingmu]=GPFit_general(D,Y,1);
[YpredDingmu, PosCovDingmu, LCLDingmu, UCLDingmu]=GPPredictX_general(Dtest,1,roptDingmu,RIDingmu,betaDingmu,RIResDingmu,sigma2Dingmu,SDDingmu,nDingmu,MinDingmu,RangeDingmu);

[roptGraepel,RIGraepel,RIResGraepel,sigma2Graepel,SDGraepel,nGraepel]=GPFit_Graepel(D,Y,1);
[YpredGraepel, PosCovGraepel, LCLGraepel, UCLGraepel]=GPPredictX_Graepel(Dtest,1,roptGraepel,RIGraepel,RIResGraepel,sigma2Graepel,SDGraepel,nGraepel);

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

res=[MAE0  MAEBmu  MAEDingmu  MAE  MAEGraepel 
RMSE0  RMSEBmu  RMSEDingmuu  RMSE RMSEGraepel 
ALPI0  ALPIBmu  ALPIDingmuu  ALPI ALPIGraepel 
Coverage0  CoverageBmu  CoverageDingmu  Coverage CoverageGraepel ]';
res
save('RESULT.mat','res')


