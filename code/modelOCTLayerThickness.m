function fVal=modelOCTLayerThickness(varargin )
% Caclulates RGC layer thickness by reference to constituent cell classes
%
% Description:
%   We aim to create a model of the thickness of the retinal ganglion cell
%   layer and how it varies with eccentricity, comparing our estimations to
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
    modelOCTLayerThickness()
%}
%{
    % Search across packing density
    modelOCTLayerThickness('forceRecalculate',true,'showPlots',false);
    myObj = @(p) modelOCTLayerThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false);
    % Set x0 to the maximum sphere packing density
    x0=0.74;
    [p,fval]=fminsearch(myObj,x0)
%}
%{
    % Search across packing density using the Drasdo midget fraction
    modelOCTLayerThickness('forceRecalculate',true,'showPlots',false);
    myObj = @(p) modelOCTLayerThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false,'midgetLinkingFuncParams',[nan]);
    % Set x0 to the maximum sphere packing density
    x0=0.74;
    [p,fval]=fminsearch(myObj,x0)
    modelOCTLayerThickness('midgetLinkingFuncParams',[nan],'packingDensity',p,'showPlots',true,'forceRecalculate',false);
%}
%{
    % Search across both packing density and midget fraction params
    modelOCTLayerThickness('forceRecalculate',true,'showPlots',false);
    myObj = @(p) modelOCTLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'objectiveType','all');
    x0=[12.0290, 0.4320   0.4171    0.9483  0.4809];
    ub=[15 0.7 0.6 1.0 0.7];
    lb=[1 0.1 0.2 0.9 0.4];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelOCTLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'objectiveType','all');
%}

%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);
p.addParameter('midgetLinkingFuncParams',[7.6048    0.4771    0.4060    0.9266],@isnumeric); % Best fit to the OCT data
p.addParameter('supportDegNasal',(0.75:.25:10)',@isnumeric);
p.addParameter('supportDegTemporal',(0.75:.25:15)',@isnumeric);
p.addParameter('packingDensity',0.4352,@isscalar);  % Best fit to the OCT data
p.addParameter('objectiveType','all',@ischar);
p.addParameter('forceRecalculate',true,@islogical);
p.addParameter('showPlots',true,@islogical);
p.addParameter('outDir','~/Desktop/KaraCloud_VSS2019_figs',@ischar);
p.addParameter('octDataFileName', ...
    fullfile(tbLocateProject('rgcPopulationModel','verbose',false),'data','rgcIplThicknessMap.mat'),...
    @ischar);

% parse
p.parse(varargin{:})


% Obtain the cell population components from the individual functions
persistent amacrine totalRGC bistratified ipRGC
if isempty(amacrine) || p.Results.forceRecalculate
    ipRGC = cell.ipRGC(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
    amacrine = cell.amacrine(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
    totalRGC = cell.totalRGC(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
    bistratified = cell.bistratified(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames, totalRGC);
end

midget = cell.midget(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames,p.Results.midgetLinkingFuncParams, totalRGC);
parasol = cell.parasol(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames, totalRGC, midget, bistratified);

%% Obtain the TOME OCT thickness measurements
persistent rgcOCTMm ratioFuncByThickness
if isempty(rgcOCTMm) || p.Results.forceRecalculate
    % Load the OCT data
    load(p.Results.octDataFileName,'rgcIplThicknessMap');
    octRadialDegreesVisualExtent = 15;
    
    % Extract the horizontal meridian
    rgciplOCTthickness = rgcIplThicknessMap(round(size(rgcIplThicknessMap,1)/2),:);
    
    % Construct a data structure to hold the OCT thickness values. Start
    % with the support in visual degrees
    rgcOCTMm.supportDeg.temporal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*octRadialDegreesVisualExtent;
    rgcOCTMm.supportDeg.nasal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*octRadialDegreesVisualExtent;
    
    % Obtain the function to express proportion of RGC+IPL thickness
    % that is RGC
    ratioFuncByThickness = rgcLayerProportion;
    
    % Obtain the total thickness of the RGC layer in mm
    rgciplOCTMm.thickMM.temporal = fliplr(rgciplOCTthickness(1:round(length(rgciplOCTthickness)/2)))./1000;
    rgciplOCTMm.thickMM.nasal = rgciplOCTthickness(round(length(rgciplOCTthickness)/2)+1:end)./1000;
    
    % Remove nans from the data and corresponding support
    notNanIdx = ~isnan(rgciplOCTMm.thickMM.temporal);
    rgciplOCTMm.thickMM.temporal = rgciplOCTMm.thickMM.temporal(notNanIdx);
    rgcOCTMm.supportDeg.temporal = rgcOCTMm.supportDeg.temporal(notNanIdx);
    
    notNanIdx = ~isnan(rgciplOCTMm.thickMM.nasal);
    rgciplOCTMm.thickMM.nasal = rgciplOCTMm.thickMM.nasal(notNanIdx);
    rgcOCTMm.supportDeg.nasal = rgcOCTMm.supportDeg.nasal(notNanIdx);
    
    % Obtain the rgc component
    rgcOCTMm.thickMM.temporal = rgciplOCTMm.thickMM.temporal .* ratioFuncByThickness(rgcOCTMm.supportDeg.temporal,rgciplOCTMm.thickMM.temporal);
    rgcOCTMm.thickMM.nasal = rgciplOCTMm.thickMM.nasal .* ratioFuncByThickness(rgcOCTMm.supportDeg.nasal,rgciplOCTMm.thickMM.nasal);
    
    % Obtain a spline fit to the thickness measurements
    for mm = [1 3]
        tmpSupport = rgcOCTMm.supportDeg.(p.Results.cardinalMeridianNames{mm})';
        tmpVals = rgcOCTMm.thickMM.(p.Results.cardinalMeridianNames{mm})';
        nonNanIdx = ~isnan(tmpVals);
        rgcOCTMm.fitDeg.(p.Results.cardinalMeridianNames{mm}) = ...
            fit( tmpSupport(nonNanIdx), tmpVals(nonNanIdx), 'smoothingspline');
    end
end


% An anonymous function for the volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

% Pick the meridian to model
cardinalMeridianNames = {'nasal','temporal'};

% Set the support. A range of 0-10 degrees matches the Curcio data.
supportDegNasal = p.Results.supportDegNasal;
supportDegTemporal = p.Results.supportDegTemporal;

for mm=1:2
    switch mm
        case 1
            supportDeg = supportDegNasal;
        case 2
            supportDeg = supportDegTemporal;
    end
    calcRGCthickness.(cardinalMeridianNames{mm}) = ((...
        amacrine.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(amacrine.diameter.fitDeg(supportDeg)) + ...
        parasol.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(parasol.diameter.fitDeg(supportDeg)) + ...
        bistratified.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(bistratified.diameter.fitDeg(supportDeg)) + ...
        midget.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(midget.diameter.fitDeg(supportDeg)) + ...
        ipRGC.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(ipRGC.diameter.fitDeg(supportDeg)) ...
        )) ./ p.Results.packingDensity ./ calc_mmSqRetina_per_degSqVisual(supportDeg,180);
end

% Objective for the nasal and temporal meridians
switch p.Results.objectiveType
    case 'all'
        nasalError = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDegNasal') - calcRGCthickness.nasal).^2));
        temporalError = sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDegTemporal') - calcRGCthickness.temporal).^2));
    case 'shape'
        nasalError = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDegNasal')./max(rgcOCTMm.fitDeg.nasal(supportDegNasal')) - calcRGCthickness.nasal./max(calcRGCthickness.nasal)).^2));
        temporalError = sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDegTemporal')./max(rgcOCTMm.fitDeg.temporal(supportDegTemporal')) - calcRGCthickness.temporal./max(calcRGCthickness.temporal)).^2));
    case 'magnitude'
        nasalError =  sqrt( (max(rgcOCTMm.fitDeg.nasal(supportDegNasal')) - max(calcRGCthickness.nasal)).^2);
        temporalError = sqrt( (max(rgcOCTMm.fitDeg.temporal(supportDegTemporal')) - max(calcRGCthickness.temporal)).^2);
end
fVal = max([nasalError temporalError]).^2;

% % Make maps of cell density
% if p.Results.makeCellMaps
%     load(p.Results.octDataFileName,'rgcIplThicknessMap');
% 
%     % Create the support for the image in degrees visual angle
%     octRadialDegreesVisualExtent = 15;
%     dim = size(rgcIplThicknessMap,1);
%     octSupportDegVis = ((1:round(dim/2))./round(dim/2)).*octRadialDegreesVisualExtent;
%     
%     % Conver the RGC+IPL thickness map to polar coordinates
%     rgcIplThicknessPolar = convertImageMapToPolarMap(rgcIplThicknessMap);
%     
%     % Obtain the function to express proportion of RGC+IPL thickness
%     % that is RGC
%     ratioFuncByThickness = rgcLayerProportion;
%     
%     
% end


%% Figure prep
if p.Results.showPlots
    
    for jj=1:4
        figHandles{jj} = figure();
    end
    
    for mm = 1:2
        meridian = cardinalMeridianNames{mm};
        
        switch mm
            case 1
                supportDeg = supportDegNasal;
                toggle = 1;
            case 2
                supportDeg = supportDegTemporal;
                toggle = -1;
        end
        
        % Plot counts
        figure(figHandles{1})
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(midget.density.fitDegSq.(meridian)(supportDeg),toggle),'-r')
        hold on
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(parasol.density.fitDegSq.(meridian)(supportDeg),toggle),'-k')
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(bistratified.density.fitDegSq.(meridian)(supportDeg),toggle),'-b')
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(amacrine.density.fitDegSq.(meridian)(supportDeg),toggle),'-g')
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(ipRGC.density.fitDegSq.(meridian)(supportDeg),toggle),'-c')
        %plot(supportDeg, midget.density.fitDegSq.(meridian)(supportDeg) + parasol.density.fitDegSq.(meridian)(supportDeg) + bistratified.density.fitDegSq.(meridian)(supportDeg) +ipRGC.density.fitDegSq.(meridian)(supportDeg) + amacrine.density.fitDegSq.(meridian)(supportDeg),'xr')
        xlabel('visual angle [deg]');
        ylabel('density [counts / sq visual deg]');
        legend({'midget','parasol','bistratified','amacrine','ipRGC'});
        xlim([-15 15]);
        pbaspect([4.58 1 1])
        box off
        grid off
        
        % Plot cell diameters
        figure(figHandles{2});
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(midget.diameter.fitDeg(supportDeg),toggle),'-r')
        hold on
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(parasol.diameter.fitDeg(supportDeg),toggle),'-k')
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(bistratified.diameter.fitDeg(supportDeg),toggle),'-b')
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(amacrine.diameter.fitDeg(supportDeg),toggle),'-g')
%        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(ipRGC.diameter.fitDeg(supportDeg),toggle),'-c')
        xlabel('visual angle [deg]');
        ylabel('soma diameter [mm]');
        legend({'midget','parasol','bistratified','amacrine'});
        xlim([-15 15]);
        pbaspect([4.58 1 1])
        box off
        grid off
        
        
        % Plot thickness
        figure(figHandles{3});
        
        plot(toggle.*toggleFlip(rgcOCTMm.supportDeg.(meridian),toggle), toggleFlip(rgcOCTMm.thickMM.(meridian),toggle), '-r');
        hold on
        plot(toggle.*toggleFlip(supportDeg,toggle), toggleFlip(calcRGCthickness.(meridian),toggle),'xk');
        legend({'RGC layer thickness OCT','Cell population model'});
        xlabel('visual angle [deg]');
        ylabel('layer thickness [mm]]');
        ylim([0 0.08]);
        xlim([-15 15]);
        pbaspect([4.58 1 1])
        box off
        grid off
    end
    
    % Plot the midget fraction model
    figure(figHandles{4});
    regularSupportPosDegVisual = 0:0.1:15;
    rgcDensitySqDegVisual = totalRGC.density.fitDegSq.temporal(regularSupportPosDegVisual');
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDegVisual );
    plot(regularSupportPosDegVisual,midgetFraction,'-r');
    hold on
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDrasdo( regularSupportPosDegVisual, rgcDensitySqDegVisual );
    plot(regularSupportPosDegVisual,midgetFraction,'--r');
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDegVisual, 'linkingFuncParams', p.Results.midgetLinkingFuncParams );
    plot(regularSupportPosDegVisual,midgetFraction,'-k');
            legend({'Dacey','Drasdo','current'});
    xlabel('visual angle [deg]')
    ylabel('Midget Fraction')
    ylim([0 1]);
    xlim([0 15]);
    pbaspect([1 2 1])

    % If the outputDir is not empty, save the figures
    if ~isempty(p.Results.outDir)
        figNames = {'modelFit_cellCounts.pdf','modelFit_cellDiameters.pdf','modelFit_cellThickness.pdf','modelFit_midgetFraction.pdf'};
        for jj=1:4
            filename = fullfile(p.Results.outDir,figNames{jj});
            print(figHandles{jj},filename,'-dpdf')
            close(figHandles{jj})
        end
    end
end

end % rgcThickness function


function output = toggleFlip(input,toggle)
if toggle==-1
    output=fliplr(input);
else
    output=input;
end

end