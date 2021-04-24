% clc,clear
% seed=sobolset(2);
% n=10;
% D_all=net(seed,2*n+2);
% Dtest=D_all(3:n+2,:);
% D=D_all(n+3:2*n+2,:);
%%
clc,clear
D=[0.0375	0.1375;
0.0625	0.4625;
0.1375	0.0125;
0.1625	0.9125;
0.2125	0.1625;
0.2625	0.6625;
0.3375	0.2625;
0.3625	0.5375;
0.4125	0.4125;
0.4625	0.8125;
0.5125	0.0625;
0.5625	0.3625;
0.6125	0.6125;
0.6875	0.9625;
0.7125	0.7375;
0.7625	0.5875;
0.8125	0.3375;
0.8625	0.8625;
0.9375	0.2125;
0.9875	0.7625];
n=size(D,1);
Duncoded=[D(:,1)*30000 D(:,2)*73+27];
Y=zeros(n,1);
for i=1:n
    Y(i)=rodexample(Duncoded(i,1),Duncoded(i,2));
end
%n=10;
Dtest=[0.0125	0.8875;
0.0875	0.6875;
0.1125	0.3125;
0.1875	0.5125;
0.2375	0.7875;
0.2875	0.3875;
0.3125	0.9375;
0.3875	0.0875;
0.4375	0.6375;
0.4875	0.2375;
0.5375	0.4875;
0.5875	0.8375;
0.6375	0.1125;
0.6625	0.2875;
0.7375	0.4375;
0.7875	0.1875;
0.8375	0.7125;
0.8875	0.0375;
0.9125	0.9875;
0.9625	0.5625];
ntest=size(Dtest,1);
Dtestuncoded=[Dtest(:,1)*30000 Dtest(:,2)*73+27];
Ytest=zeros(ntest,1);
for i=1:ntest
    Ytest(i)=rodexample(Dtestuncoded(i,1),Dtestuncoded(i,2));
end
%%
[ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range]=GPFit0(D,Y,1);
[Ypred0, PosCov0, LCL0, UCL0]=GPPredictX0(Dtest,1,ropt0,RI0,beta0,RIRes0,sigma20,SD0,n0,Min,Range);

[alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k]=BMGPFit(D,Y,1);
[Ypred, PosCov, LCL, UCL]=BMGPPredictX(Dtest,1,alphaopt,etaopt,ropt,a0opt,s2,unscaledsigmaD,QI,QIRes,Min,Range,SD,n,k);

[alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k]=BMGPFit_general(D,Y,0);
[YpredBmu, PosCovBmu, LCLBmu, UCLBmu]=BMGPPredictX_general(Dtest,0,alphaoptBmu,roptBmu,a0optBmu,s2Bmu,QIBmu,QIResBmu,SDBmu,n,k);

[alphaoptDingmu,roptDingmu,a0optDingmu,s2Dingmu,QIDingmu,QIResDingmu,SDDingmu,n,k]=BMGPFit_general(D,Y,1);
[YpredDingmu, PosCovDingmu, LCLDingmu, UCLDingmu]=BMGPPredictX_general(Dtest,1,alphaoptDingmu,roptDingmu,a0optDingmu,s2Dingmu,QIDingmu,QIResDingmu,SDDingmu,n,k);

%%
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

res=[MAE0  MAEBmu  MAEDingmu  MAE 
RMSE0  RMSEBmu  RMSEDingmuu  RMSE
ALPI0  ALPIBmu  ALPIDingmuu  ALPI
Coverage0  CoverageBmu  CoverageDingmu  Coverage]';
save('RESULT.mat','res')


