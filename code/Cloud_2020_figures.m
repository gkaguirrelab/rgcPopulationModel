
%% Search across packing density
% This is a search using the default Dacey midget fraction model
args = {};
modelGCLayerThickness('showPlots',false,'forceRecalculate',true, args{:});
myObj = @(p) modelGCLayerThickness('packingDensity',p(1),'showPlots',false,'forceRecalculate',false, args{:});
x0=[0.6122];
ub=[0.9];
lb=[0.1];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('packingDensity',p(1),'showPlots',true,'forceRecalculate',true, args{:});


% This is a search using the default Dacey midget fraction model and Liu
% cell densities
args = {'totalRGCSource','Liu','meridianSetName','horiz'};
modelGCLayerThickness('showPlots',false,'forceRecalculate',true, args{:});
myObj = @(p) modelGCLayerThickness('packingDensity',p(1),'showPlots',false,'forceRecalculate',false, args{:});
x0=[0.6122];
ub=[0.9];
lb=[0.1];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('packingDensity',p(1),'showPlots',true,'forceRecalculate',true, args{:});



% This is a search using the Drasdo midget fraction model
modelGCLayerThickness('midgetLinkingFuncParams',[1.0000   41.0300    0.0000    0.8928],'showPlots',false,'forceRecalculate',true);
myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',[1.0000   41.0300    0.0000    0.8928],'packingDensity',p(1),'showPlots',false,'forceRecalculate',false);
x0=[0.6122];
ub=[0.9];
lb=[0.1];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('midgetLinkingFuncParams',[1.0000   41.0300    0.0000    0.8928],'packingDensity',p(1),'showPlots',true,'forceRecalculate',true);


% Search over full model
args = {'midgetModel','fixed'};
modelGCLayerThickness('showPlots',false,'forceRecalculate',true,args{:});
myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'cellSizeParams',p(6:7), 'showPlots',false,'forceRecalculate',false,args{:});
x0=[6.8151   21.4511    0.4668    0.9577  0.5336 0, 0];
ub=[10, 50, 0.5, 1.000, 0.8, 0, 0];
lb=[1, 10, 0.4, 0.9, 0.1, 0, 0];
[p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'cellSizeParams',p(6:7),'showPlots',true,'forceRecalculate',true,args{:});



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
totalRGC = cell.totalRGC_Curcio;
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

