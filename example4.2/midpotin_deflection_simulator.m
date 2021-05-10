function y=midpotin_deflection_simulator(x)
numberOfPDE = 2;
model = createpde(numberOfPDE);
E = 1.0e6; % modulus of elasticity
nu = .3; % Poisson's ratio
len = 10.0; % side length for the square plate
D = E*x^3/(12*(1 - nu^2));
hmax = len/20; % mesh size parameter
pres = 2; % external pressure
gdm = [3 4 0 len len 0 0 0 len len]';
g = decsg(gdm,'S1',('S1')');
geometryFromEdges(model,g);
c = [1 0 1 D 0 D]';
a = [0 0 1 0]';
f = [0 pres]';
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);
k = 1e7;
bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:4),...
'g',[0 0],'q',[0 0; k 0]);
generateMesh(model, 'Hmax', hmax);
res = solvepde(model);
u = res.NodalSolution;
numNodes = size(model.Mesh.Nodes,2);
y=interpolateSolution(res,len/2,len/2,1);
end

