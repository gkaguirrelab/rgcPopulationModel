function [fVal, figHandles] = modelGCLayerThickness(varargin )
% Caclulates GC layer thickness by reference to constituent cell classes
%
% Description:
%   We aim to create a model of the thickness of the retinal ganglion cell
%   layer and how it varies with eccentricity, comparing our estimates to
%   reports from histology and OCT scanned images. We will calculate the
%   thickness of the RGC layer using data from histology on the population
%   of midget, parasol, and bistratified ganglion cells as well as
%   displaced amacrine cells. Specifically, the proportion of cell types,
%   their size, and their densities and how these three amounts vary with
%   eccentricity. The volume of space each cell type takes up will be
%   combined to create a coherent digital map of the human retinal ganglion
%   cell layer. This project utilizes largely generalized averages of cell
%   sizes and proportions, and does not take into account the different
%   subtypes of ganglion cells such as ON and OFF center midget and parasol
%   cells, widefield, smallfield, and diffuse bistratified cells, or the
%   many types of amacrine cells. The use of average proportions, sizes,
%   and densities, are noted throughout and still create an accurate
%   estimate of RGC thickness.
%
% Inputs:
%
% Optional key / value pairs:
%  'midgetLinkingFuncParams' - 1x4 vector. The parameters that define the
%                           form of the midget fraction function.
%  'packingDensity'       - Scalar. We model the cells as spheres. For a
%                           single sphere size, the most dense possible
%                           packing is Keppler's limit of ~0.74. In the RGC
%                           there are multiple cell classes. Except in the
%                           far periphery the two most numerous types of
%                           cells are the midget and parasol. If smaller
%                           cells are more numerous than larger cells, it
%                           is possible in principle to pack small cells
%                           into the spaces between larger cells, and
%                           achieve a density higher than Kepler's limit.
%                           When the ratio of the smaller to larger cell is
%                           greater than ~0.41, then the smaller is too
%                           large to pack into the larger. This is the
%                           situation we face in this model, so we assume
%                           Keppler's limit for sphere packing.
%  'cellSizeParams'       - 1x2 vector. While there is some evidence that
%                           cell size increases with eccentricity, there
%                           are limited data to quantify this property. We
%                           model soma size for each class of RGC by
%                           starting with the mean of values reported
%                           across a range of ecentricities. This size
%                           varies in proportion linearly across
%                           eccentricity under the control of two
%                           parameters which are in common for all cell
%                           classes.
%  'forceRecalculate'     - When set to true forces the routine to
%                           re-calculate the thickness ratio function from
%                           the source data.
%  'midgetModel'          - Char vector. Valid values are
%                               {'fixed','proportional'}
%                           which switches if the midget fraction varies
%                           simply as a function of eccentricity, or
%                           following a more complicated model of
%                           proportion of total RGCs to that point.
%  'totalRGCSource'       - Char vector. Valid values are
%                               {'Curcio','Liu'}
%                           which determines the source of the values of
%                           total RGCs along the meridians.
%  'meridianSetName'      - Char vector. Valid values are:
%                               {'horiz','vert','both'}
%                           which determines the set of cardinal meridians
%                           along which the calculations are performed.
%                           For the 'Liu' totalRGCSource, only the horiz
%                           set is available.
%  'forceRecalculate'     - Force the routine to re-calculate the cell
%                           populations. Can be set to false to speed
%                           repeated calls to the function.
%
% Outputs:
%
% Examples:
%{
    % Simple example
    modelGCLayerThickness()
%}
%{
    % Search across packing density
    myObj = @(p) modelGCLayerThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false);
    x0=[0.6122];
    ub=[0.9];
    lb=[0.1];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelGCLayerThickness('packingDensity',p,'showPlots',true,'forceRecalculate',false);
%}


%% input parser
p = inputParser;

% Optional analysis params
%p.addParameter('midgetLinkingFuncParams',[ 12.0890    0.4335   0.5 0.95],@isnumeric); % Values for the midget proportional model
p.addParameter('midgetLinkingFuncParams',[ 6.8151   21.4511    0.4668    0.9577 ],@isvector); % Values for the midget fixed model
p.addParameter('packingDensity',0.6214,@isscalar);  % Best fit to the OCT data
p.addParameter('cellSizeParams',[0,0],@isvector);  % Best fit to the OCT data
p.addParameter('midgetModel','fixed',@ischar); % Valid values are {'fixed','proportional'}
p.addParameter('totalRGCSource','Curcio',@ischar); % Valid values are {'Curcio','Liu'}
p.addParameter('totalRGCSupportShift',-0.33,@isscalar); % Valid values are {'Curcio','Liu'}
p.addParameter('meridianSetName','both',@ischar);
p.addParameter('forceRecalculate',true,@islogical);
p.addParameter('showPlots',true,@islogical);
p.addParameter('showInputPlots',false,@islogical);
p.addParameter('outDir','~/Desktop/KaraCloud_VSS2019_figs',@ischar);

% parse
p.parse(varargin{:})

% Obtain the ganglion cell thickness data
thickData = prepareGCThickData(p.Results.meridianSetName);

% Obtain the mean ganglion cell size data
sizeData = prepareMeanCellSizeData_Liu(p.Results.meridianSetName);

% Obtain the cell population components from the individual functions
% This first set are invariant with respect to the midget fraction and the
% cellSizeParams
persistent totalRGC ipRGC
if isempty(totalRGC) || p.Results.forceRecalculate
    ipRGC = cell.ipRGC(p.Results.showInputPlots);
    switch p.Results.totalRGCSource
        case 'Curcio'
            totalRGC = cell.totalRGC_Curcio(p.Results.totalRGCSupportShift,p.Results.showInputPlots);
        case 'Liu'
            totalRGC = cell.totalRGC_Liu(p.Results.showInputPlots);
    end
end

% These vary with the midget fraction linking parameters
switch p.Results.midgetModel
    case 'proportional'
        midget = cell.midget_proportional(totalRGC,p.Results.midgetLinkingFuncParams,p.Results.cellSizeParams,p.Results.showInputPlots);
    case 'fixed'
        midget = cell.midget_fixed(totalRGC,p.Results.midgetLinkingFuncParams,p.Results.cellSizeParams,p.Results.showInputPlots);
end

% These vary with the cellSizeParams, so must be re-calculated
amacrine = cell.amacrine(p.Results.cellSizeParams,p.Results.showInputPlots);
bistratified = cell.bistratified(totalRGC,p.Results.cellSizeParams,p.Results.showInputPlots);
parasol = cell.parasol(totalRGC, midget,bistratified,p.Results.cellSizeParams, p.Results.showInputPlots);

% An anonymous function for the volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

% Loop over the meridians provided in the data. Build the model and obtain
% the error.
thickError = [];
sizeError = [];

for mm=1:length(thickData)
    
    % Confirm that the label matches for the thick and size data
    if ~strcmp(thickData(mm).label,sizeData(mm).label)
        error('The data arrays are not index aligned')
    end
    
    % Prepare this entry of the models
    thickModel(mm).label = thickData(mm).label;
    thickModel(mm).angle = thickData(mm).angle;
    
    sizeModel(mm).label = sizeData(mm).label;
    sizeModel(mm).angle = sizeData(mm).angle;
    
    % Find the matching cell component entry
    idx = strcmp([totalRGC(:).label],thickData(mm).label);
    
    % The model thickness is the sum of the layer components
    thickModel(mm).thickMM = @(x) ((...
        amacrine(idx).countsDegSq(x) .* sVol(amacrine(idx).diameter(x)) + ...
        parasol(idx).countsDegSq(x) .* sVol(parasol(idx).diameter(x)) + ...
        bistratified(idx).countsDegSq(x) .* sVol(bistratified(idx).diameter(x)) + ...
        midget(idx).countsDegSq(x) .* sVol(midget(idx).diameter(x)) + ...
        ipRGC(idx).countsDegSq(x) .* sVol(ipRGC(idx).diameter(x)) ...
        )) ./ p.Results.packingDensity ./ calc_mmSqRetina_per_degSqVisual(x',totalRGC(mm).angle);
    
    % The modeled mean cell size is the weighted mean of the cell sizes
    sizeModel(mm).diameter = @(x)  ((...
        amacrine(idx).countsDegSq(x) .* amacrine(idx).diameter(x) + ...
        parasol(idx).countsDegSq(x) .* parasol(idx).diameter(x) + ...
        bistratified(idx).countsDegSq(x) .* bistratified(idx).diameter(x) + ...
        midget(idx).countsDegSq(x) .* midget(idx).diameter(x) + ...
        ipRGC(idx).countsDegSq(x) .* ipRGC(idx).diameter(x) ...
        )) ./ ((...
        amacrine(idx).countsDegSq(x) + ...
        parasol(idx).countsDegSq(x) + ...
        bistratified(idx).countsDegSq(x) + ...
        midget(idx).countsDegSq(x) + ...
        ipRGC(idx).countsDegSq(x) ...
        ));
        
    % Evaluate the model at the data locations
    dataModelThickDif = thickData(mm).thickMM - thickModel(mm).thickMM(thickData(mm).supportDeg);
    dataModelSizeDif = sizeData(mm).diameter - sizeModel(mm).diameter(sizeData(mm).supportDeg);
    
    thickError = [thickError; dataModelThickDif];
    sizeError = [sizeError; dataModelSizeDif];
    
end

% Calculate the error metric as RMSE in fiting thickness
fVal = sqrt(nansum(thickError.^2)) .* (1+sqrt(nansum(sizeError.^2)));


%% Show some plots 
if p.Results.showPlots
    
    % Plot the thickness model fit
    figHandles(1) = figure();
    modelSupportDeg = 0:0.1:25;
    for mm = 1:length(thickData)
        xVals = thickData(mm).supportDeg;
        yVals = thickData(mm).thickMM;
        plotColor = [0.75 0.75 0.75];
        switch thickData(mm).label
            case 'temporal'
                subplot(2,1,1)
                xVals = -xVals;
                xValsModel = -modelSupportDeg;
            case 'nasal'
                subplot(2,1,1)
                xValsModel = modelSupportDeg;
            case 'inferior'
                subplot(2,1,2)
                xVals = -xVals;
                xValsModel = -modelSupportDeg;
            case 'superior'
                subplot(2,1,2)
                xValsModel = modelSupportDeg;
        end
        area(xVals,yVals,'FaceColor',plotColor,'EdgeColor','none')
        hold on
        xlim([-25 25]);
        ylim([0 0.075]);
        plot(xValsModel,thickModel(mm).thickMM(modelSupportDeg),'-','Color',[211 60 33]./255,'LineWidth',3);
        title(['GC thickness - ' thickModel(mm).label])
        drawnow
    end
    drawnow
    
    % Plot the mean size model fit
    figHandles(2) = figure();
    supportDeg = 1:0.1:40;
    for mm = 1:length(sizeData)
        subplot(2,2,mm)
        plot(sizeData(mm).supportDeg,sizeData(mm).diameter,'xk');
        hold on
        plot(supportDeg,sizeModel(mm).diameter(supportDeg),'-r');
        title(['Mean GC size - ' sizeModel(mm).label])
        ylim([0 0.02]);
        xlim([0 40]);
        drawnow
    end
    
    % Plot the midget fraction model
    figHandles(3) = figure();
    supportDeg = 0:0.1:40;
    idx = strcmp([totalRGC(:).label],'temporal');
    countsDegSqTotal = totalRGC(idx).countsDegSq(supportDeg);
    [~, daceyMidgetFraction, daceyDataSupportPosDegVisual] = calcDaceyMidgetFractionByEccenDegVisual(supportDeg);
    plot(daceyDataSupportPosDegVisual,daceyMidgetFraction,'or');
    hold on
    drasdoMidgetFraction = calcDrasdoMidgetFractionByVisualEccen(supportDeg,0.8928,41.03);
    plot(supportDeg,drasdoMidgetFraction,'--r');
    
    switch p.Results.midgetModel
        case 'proportional'
            [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( supportDeg, countsDegSqTotal,'linkingFuncParams',p.Results.midgetLinkingFuncParams(1:2),'minMidgetFractionRatio',p.Results.midgetLinkingFuncParams(3),'maxMidgetFractionRatio',p.Results.midgetLinkingFuncParams(4));
        case 'fixed'
            [~, logisticFunc] = daceyMidgetFractionByEccenDegVisual();
            midgetFraction = logisticFunc(...
                p.Results.midgetLinkingFuncParams(1),...
                p.Results.midgetLinkingFuncParams(2),...
                p.Results.midgetLinkingFuncParams(3),...
                p.Results.midgetLinkingFuncParams(4),...
                supportDeg);
    end
    plot(supportDeg,midgetFraction,'-k');
    legend({'Dacey','Drasdo','current'},'Location','southwest');
    xlabel('visual angle [deg]')
    ylabel('Midget Fraction')
    ylim([0.4 1]);
    xlim([0 40]);
    pbaspect([1 2 1])
    
    % Plot the cell densities
    supportDeg = 0:0.1:25;
    figHandles(4) = figure();
    cellTypes = {'midget','parasol','bistratified','amacrine'};
    plotColors = {'-r','-k','-b','-g'};
    for cc = 1:length(cellTypes)
        statement = ['yVals = ' cellTypes{cc} '(idx).countsDegSq(supportDeg);'];
        eval(statement);
        plot(supportDeg,yVals,plotColors{cc});
        hold on
    end
    ylim([0 2000]);
    ylabel('Cell density [counts/deg^2]');
    xlabel('Eccentricity [deg]')
end

end

