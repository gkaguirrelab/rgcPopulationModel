function fVal=modelGCLayerThickness(varargin )
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
%  'forceRecalculate'     - When set to true forces the routine to
%                           re-calculate the thickness ratio function from
%                           the source data.

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
%{
    % Search across model params to fit the data
    % FIXED midget model
    myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'midgetModel','fixed');
    x0=[6.2665,21.9499,0.5,0.9500, 0.6];
    ub=[10, 30, 0.6, 1.000, 1.0];
    lb=[1, 10, 0.4, 0.925, 0.4];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'midgetModel','fixed');
%}
%{
    % Search across model params to fit the data
    % PROPORTIONAL midget model
    myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'midgetModel','proportional');
    x0=[12, 0.4, 0.5, 0.950, 0.6];
    ub=[30, 0.6, 0.6, 1.000, 1.0];
    lb=[05, 0.2, 0.4, 0.925, 0.4];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'midgetModel','proportional');
%}
%{
    % Search across model params to fit the data
    % FIXED midget model
    myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'midgetModel','fixed');
    x0=[6.2665,21.9499,0.5,0.95, 0.6];
    ub=[10, 30, 0.5, 1.000, 1.0];
    lb=[1, 10, 0.4, 0.9, 0.4];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'midgetModel','fixed');
%}


%% input parser
p = inputParser;

% Optional analysis params
%p.addParameter('midgetLinkingFuncParams',[ 12.0890    0.4335   0.5 0.95],@isnumeric); % Values for the midget proportional model
p.addParameter('midgetLinkingFuncParams',[ 6.2665,21.9499,0.45,0.96],@isnumeric); % Values for the midget fixed model
p.addParameter('packingDensity',0.6214,@isscalar);  % Best fit to the OCT data
p.addParameter('forceRecalculate',true,@islogical);
p.addParameter('showPlots',true,@islogical);
p.addParameter('meridianSetName','both',@ischar);
p.addParameter('midgetModel','fixed',@ischar); % Valid values are {'fixed','proportional'}
p.addParameter('outDir','~/Desktop/KaraCloud_VSS2019_figs',@ischar);

% parse
p.parse(varargin{:})

% Obtain the ganglion cell thickness data
thickData = prepareGCThickData(p.Results.meridianSetName);

% Obtain the mean ganglion cell size data
sizeData = prepareMeanCellSizeData(p.Results.meridianSetName);

% Obtain the cell population components from the individual functions
% This first set are invariant with respect to the midget fraction
persistent amacrine totalRGC bistratified ipRGC
if isempty(amacrine) || p.Results.forceRecalculate
    ipRGC = cell.ipRGC();
    amacrine = cell.amacrine();
    totalRGC = cell.totalRGC();
    bistratified = cell.bistratified(totalRGC);
end

% These vary with the midget fraction linking parameters
switch p.Results.midgetModel
    case 'proportional'
        midget = cell.midget_proportional(totalRGC, p.Results.midgetLinkingFuncParams);
    case 'fixed'
        midget = cell.midget_fixed(totalRGC, p.Results.midgetLinkingFuncParams);
end
parasol = cell.parasol(totalRGC, midget, bistratified);

% An anonymous function for the volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

% Loop over the meridians provided in the data. Build the model and obtain
% the error.
fValThick = 0;
fValSize = 0;
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
    
    fValThick = fValThick + sqrt(nansum(dataModelThickDif.^2));
    fValSize = fValSize + sqrt(nansum(dataModelSizeDif.^2));
    
    fVal = fValThick; % + fValSize/10;
end

if p.Results.showPlots
    
    % Plot the thickness model fit
    figure
    supportDeg = 0:0.1:25;
    for mm = 1:length(thickData)
        subplot(2,2,mm)
        plot(thickData(mm).supportDeg,thickData(mm).thickMM,'xk');
        hold on
        plot(supportDeg,thickModel(mm).thickMM(supportDeg),'-r');
        title(['GC thickness - ' thickModel(mm).label])
        drawnow
    end
    
    % Plot the mean size model fit
    figure
    supportDeg = 0:0.1:40;
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
    figure
    supportDeg = 0:0.1:40;
    countsDegSqTotal = totalRGC(3).countsDegSq(supportDeg);
    [~, daceyMidgetFraction, daceyDataSupportPosDegVisual] = calcDaceyMidgetFractionByEccenDegVisual(supportDeg);
    plot(daceyDataSupportPosDegVisual,daceyMidgetFraction,'or');
    hold on
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDrasdo( supportDeg, countsDegSqTotal );
    plot(supportDeg,midgetFraction,'--r');

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
end

end

