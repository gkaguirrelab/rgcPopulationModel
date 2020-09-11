
%% Search across packing density
% This is a search using the default Dacey midget fraction model
myObj = @(p) modelGCLayerThickness('packingDensity',p(1),'showPlots',false,'forceRecalculate',false);
x0=[0.6122];
ub=[0.9];
lb=[0.1];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('packingDensity',p(1),'showPlots',true,'forceRecalculate',true);

% This is a search using the Drasdo midget fraction model
myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',[1.0000   41.0300    0.0000    0.8928],'packingDensity',p(1),'showPlots',false,'forceRecalculate',false);
x0=[0.6122];
ub=[0.9];
lb=[0.1];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('midgetLinkingFuncParams',[1.0000   41.0300    0.0000    0.8928],'packingDensity',p(1),'showPlots',true,'forceRecalculate',true);


% Starting with the Drasdo x0 midget fraction, the model finds something
% close to Dacey
myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'cellSizeSlope',p(6), 'showPlots',false,'forceRecalculate',false,'midgetModel','fixed');
x0=[1.0000   41.0300    0.0000    0.8928, 0.5336, 0];
ub=[10, 50, 0.5, 1.000, 1.0, 0];
lb=[1, 10, 0, 0.7, 0.1, 0];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'cellSizeSlope',p(6),'showPlots',true,'forceRecalculate',false,'midgetModel','fixed');

% Create a plot of GCL thickness
figure
thickData = prepareGCThickData('both');
for ii = 1:4
    xVals = thickData(ii).supportDeg;
    yVals = thickData(ii).thickMM;
    plotColor = [0.75 0.75 0.75];
    switch thickData(ii).label
        case 'temporal'
            subplot(2,1,1)
            xVals = -xVals;
        case 'nasal'
            subplot(2,1,1)
        case 'inferior'
            subplot(2,1,2)
            xVals = -xVals;
        case 'superior'
            subplot(2,1,2)
    end
    area(xVals,yVals,'FaceColor',plotColor,'EdgeColor','none')
    xlim([-25 25]);
    ylim([0 0.075]);
    hold on
end


% Create a plot of total GC counts
totalRGC = cell.totalRGC;
for ii = 1:4
    xVals = 0:0.1:25;
    yVals = totalRGC(ii).countsDegSq(xVals);
    plotColor = [211 60 33]./255;
    switch totalRGC(ii).label{1}
        case 'temporal'
            subplot(2,1,1)
yyaxis right
            xVals = -xVals;
        case 'nasal'
            subplot(2,1,1)
yyaxis right
        case 'inferior'
            subplot(2,1,2)
yyaxis right
            xVals = -xVals;
        case 'superior'
            subplot(2,1,2)
yyaxis right
    end
    plot(xVals,yVals,'-','Color',plotColor,'LineWidth',3)
    xlim([-25 25]);
    ylim([0 2750]);
    hold on
end

