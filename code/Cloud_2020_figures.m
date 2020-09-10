
%% Search across packing density
% This is a search using the default Dacey midget fraction model
myObj = @(p) modelGCLayerThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false);
x0=0.6122;
ub=0.9;
lb=0.1;
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub);
p = 0.5303;
modelGCLayerThickness('packingDensity',p,'showPlots',true,'forceRecalculate',false);


myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'midgetModel','fixed');
x0=[6.2665,21.9499,0.5,0.95, 0.6];
ub=[10, 30, 0.5, 1.000, 1.0];
lb=[1, 10, 0.4, 0.9, 0.4];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub);
p = [9.9884   26.6404    0.4945    0.9000    0.6179];
modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'midgetModel','fixed');
