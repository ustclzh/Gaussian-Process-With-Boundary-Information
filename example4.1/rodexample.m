function temp=rodexample(R,Tleft)
%R in (0,50000] 
%Tleft in [27,100]
thermalModelT = createpde('thermal','transient');
g = decsg([3 4 -1.5 1.5 1.5 -1.5 0 0 .2 .2]');
geometryFromEdges(thermalModelT,g);
%pdegplot(thermalModelT,'EdgeLabels','on');
generateMesh(thermalModelT,'Hmax',0.1);
rho = 7800; % density, kg/m^3
cp = 500; % specific heat, W-s/(kg-degree C)
k = rho*cp/R; % thermal reffusivity, s/m^2
q = 0; % heat source, W/m^3
thermalProperties(thermalModelT,'ThermalConductivity',k,'MassDensity',rho,'SpecificHeat',cp);
internalHeatSource(thermalModelT,q);
thermalBC(thermalModelT,'Edge',2,'HeatFlux',0);
thermalBC(thermalModelT,'Edge',[1 3],'HeatFlux',0);
thermalBC(thermalModelT,'Edge',4,'Temperature',Tleft);
tfinal = 86400;
tlist = 0:60:tfinal;
thermalIC(thermalModelT,27);
thermalModelT.SolverOptions.ReportStatistics = 'on';
result = solve(thermalModelT,tlist);
T = result.Temperature;
%figure;
%pdeplot(thermalModelT,'XYData',T(:,end),'Contour','on');
%axis equal
%title(sprintf('Transient Temperature at Final Time (%g seconds)',tfinal));
p = thermalModelT.Mesh.Nodes;
nodesRightEnd = find(p(1,:) > 1.5-eps);
nodeCenter = nodesRightEnd(p(2,nodesRightEnd) < eps)
%nodeOuter = nodesRightEnd(p(2,nodesRightEnd) > 0.2-eps)
%figure;
%plot(tlist,T(nodeCenter,:));
%hold all
%plot(tlist,T(nodeOuter,:),'--');
%title 'Temperature at Left End as a Function of Time'
%xlabel 'Time, seconds'
%ylabel 'Temperature, degrees-C'
%grid on;
%legend('Center Axis','Outer Surface','Location','SouthEast');
temp=T(nodeCenter,end);