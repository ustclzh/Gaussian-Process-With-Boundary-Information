function temp=rodexample(h,Tleft)
%%
global halfl radius
%h in [0,50] 
%Tleft in [20,100]

thermalModelT = createpde('thermal','steadystate');
g = decsg([3 4 -halfl halfl halfl -halfl 0 0 radius radius]');
geometryFromEdges(thermalModelT,g);
%pdegplot(thermalModelT,'EdgeLabels','on');
generateMesh(thermalModelT,'Hmax',min(halfl/10,radius/2));
q = 0; % heat source, W/m^3
kFunc = @(region,state) (31033./(184.86+state.u)).*region.y; 
%https://www.efunda.com/materials/elements/TC_Table.cfm?Element_ID=Si, 
%rational polynomial fit to data from 250K to 1000K
qFunc = @(region,state) q*region.y;
thermalProperties(thermalModelT,'ThermalConductivity',kFunc);
internalHeatSource(thermalModelT,qFunc);
thermalBC(thermalModelT,'Edge',2,'HeatFlux',0);
outerCC = @(region,~) h*region.y;
thermalBC(thermalModelT,'Edge',3,'ConvectionCoefficient',outerCC,'AmbientTemperature',20);
thermalBC(thermalModelT,'Edge',4,'Temperature',Tleft);
thermalModelT.SolverOptions.ReportStatistics = 'off';
result = solve(thermalModelT);
T = result.Temperature;
p = thermalModelT.Mesh.Nodes;
nodesRightEnd = find(p(1,:) > halfl-eps);
nodeCenter = nodesRightEnd(p(2,nodesRightEnd) < eps);
temp=T(nodeCenter,end);