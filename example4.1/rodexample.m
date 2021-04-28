function temp=rodexample(h,Tleft)
%%
global halfl radius
%h in [0,50] 
%Tleft in [27,100]

thermalModelT = createpde('thermal','steadystate');
g = decsg([3 4 -halfl halfl halfl -halfl 0 0 radius radius]');
geometryFromEdges(thermalModelT,g);
%pdegplot(thermalModelT,'EdgeLabels','on');
generateMesh(thermalModelT,'Hmax',0.1);
%rho = 7800; % density, kg/m^3
%cp = 500; % specific heat, W-s/(kg-degree C)
q = 0; % heat source, W/m^3
kFunc = @(region,state) (54-3.33*10^-2*state.u).*region.y;
qFunc = @(region,state) q*region.y;
thermalProperties(thermalModelT,'ThermalConductivity',kFunc);
internalHeatSource(thermalModelT,qFunc);
thermalBC(thermalModelT,'Edge',2,'HeatFlux',0);
outerCC = @(region,~) h*region.y;
thermalBC(thermalModelT,'Edge',3,'ConvectionCoefficient',outerCC,'AmbientTemperature',27);
thermalBC(thermalModelT,'Edge',4,'Temperature',Tleft);
thermalModelT.SolverOptions.ReportStatistics = 'off';
result = solve(thermalModelT);
T = result.Temperature;
p = thermalModelT.Mesh.Nodes;
nodesRightEnd = find(p(1,:) > halfl-eps);
nodeCenter = nodesRightEnd(p(2,nodesRightEnd) < eps);
temp=T(nodeCenter,end);