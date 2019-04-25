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
    modelOCTLayerThickness('makeCellMaps',false)
%}
%{
    % Search across midget fraction linking params
    myObj = @(p) modelOCTLayerThickness('midgetLinkingFuncParams',p,'showPlots',false,'forceRecalculate',false,'objectiveType','shape','supportDeg',(1:1:10)');
    % Set x0 to be the parameters for the Dacey data
    x0=[12.0290 0.4320];
    ub=[14 0.6];
    lb=[10 0.2];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelOCTLayerThickness('midgetLinkingFuncParams',p,'showPlots',true,'forceRecalculate',false,'objectiveType','shape');
%}
%{
    % Search across packing density
    myObj = @(p) modelOCTLayerThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false,'objectiveType','magnitude','supportDeg',(1:1:10)');
    % Set x0 to the maximum sphere packing density
    x0=0.74;
    [p,fval]=fminsearch(myObj,x0)
%}
%{
    % Search across both packing density and midget fraction params
    myObj = @(p) modelOCTLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false,'objectiveType','all','supportDeg',(1:1:10)');
    x0=[12.0290 0.4320  0.41  0.95  0.5556];
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
p.addParameter('midgetLinkingFuncParams',[7.1901    0.4269    0.4072    0.9443],@isnumeric); % Best fit to the OCT data
p.addParameter('supportDeg',(0:1:10)',@isnumeric);
p.addParameter('packingDensity',0.5858,@isscalar);  % Best fit to the OCT data
p.addParameter('objectiveType','all',@ischar);
p.addParameter('forceRecalculate',true,@islogical);
p.addParameter('showPlots',true,@islogical);
p.addParameter('makeCellMaps',false,@islogical);
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
supportDeg = p.Results.supportDeg;

for mm=1:2
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
        nasalError = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDeg') - calcRGCthickness.nasal).^2));
        temporalError = sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDeg') - calcRGCthickness.temporal).^2));
    case 'shape'
        nasalError = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDeg')./max(rgcOCTMm.fitDeg.nasal(supportDeg')) - calcRGCthickness.nasal./max(calcRGCthickness.nasal)).^2));
        temporalError = sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDeg')./max(rgcOCTMm.fitDeg.temporal(supportDeg')) - calcRGCthickness.temporal./max(calcRGCthickness.temporal)).^2));
    case 'magnitude'
        nasalError =  sqrt( (max(rgcOCTMm.fitDeg.nasal(supportDeg')) - max(calcRGCthickness.nasal)).^2);
        temporalError = sqrt( (max(rgcOCTMm.fitDeg.temporal(supportDeg')) - max(calcRGCthickness.temporal)).^2);
end
fVal = max([nasalError temporalError]).^2;

% Make maps of cell density
if p.Results.makeCellMaps
    load(p.Results.octDataFileName,'rgcIplThicknessMap');

    % Create the support for the image in degrees visual angle
    octRadialDegreesVisualExtent = 15;
    dim = size(rgcIplThicknessMap,1);
    octSupportDegVis = ((1:round(dim/2))./round(dim/2)).*octRadialDegreesVisualExtent;
    
    % Conver the RGC+IPL thickness map to polar coordinates
    rgcIplThicknessPolar = convertImageMapToPolarMap(rgcIplThicknessMap);
    
    % Obtain the function to express proportion of RGC+IPL thickness
    % that is RGC
    ratioFuncByThickness = rgcLayerProportion;
    
    
end


%% Figure prep
if p.Results.showPlots
    for mm = 1:2
        figure
        meridian = cardinalMeridianNames{mm};
        
        % Plot counts
        subplot(2,2,1)
        plot(supportDeg, midget.density.fitDegSq.(meridian)(supportDeg))
        hold on
        plot(supportDeg, parasol.density.fitDegSq.(meridian)(supportDeg))
        plot(supportDeg, bistratified.density.fitDegSq.(meridian)(supportDeg))
        plot(supportDeg, amacrine.density.fitDegSq.(meridian)(supportDeg))
        plot(supportDeg, ipRGC.density.fitDegSq.(meridian)(supportDeg))
        %plot(supportDeg, midget.density.fitDegSq.(meridian)(supportDeg) + parasol.density.fitDegSq.(meridian)(supportDeg) + bistratified.density.fitDegSq.(meridian)(supportDeg) +ipRGC.density.fitDegSq.(meridian)(supportDeg) + amacrine.density.fitDegSq.(meridian)(supportDeg),'xr')
        xlabel('eccentricity [deg visual]');
        ylabel('density [counts / sq visual deg]');
        legend({'midget','parasol','bistratified','amacrine','ipRGC'});
        
        % Plot cell volumes
        subplot(2,2,2)
        plot(supportDeg, sVol(midget.diameter.fitDeg(supportDeg)))
        hold on
        plot(supportDeg, sVol(parasol.diameter.fitDeg(supportDeg)))
        plot(supportDeg, sVol(bistratified.diameter.fitDeg(supportDeg)))
        plot(supportDeg, sVol(amacrine.diameter.fitDeg(supportDeg)))
        plot(supportDeg, sVol(ipRGC.diameter.fitDeg(supportDeg)))
        xlabel('eccentricity [deg visual]');
        ylabel('individual cell volume [mm^3]');
        legend({'midget','parasol','bistratified','amacrine','ipRGC'});
        
        
        % Plot thickness
        subplot(2,2,3)
        
        plot(rgcOCTMm.supportDeg.(meridian), rgcOCTMm.thickMM.(meridian), 'xk');
        hold on
        plot(supportDeg, calcRGCthickness.(meridian),'-k');
        legend({'RGC layer thickness OCT','Cell population model'});
        xlabel('eccentricity [deg visual]');
        ylabel('layer thickness [mm]]');
        ylim([0 0.08]);
        xlim([0 15]);
    end
    
    % Plot the midget fraction model
    regularSupportPosDegVisual = 0:0.1:30;
    rgcDensitySqDegVisual = totalRGC.density.fitDegSq.temporal(regularSupportPosDegVisual');
    figure
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDegVisual );
    plot(regularSupportPosDegVisual,midgetFraction,'-k');
    hold on
    [ ~, midgetFraction ] = transformRGCToMidgetRGCDensityDacey( regularSupportPosDegVisual, rgcDensitySqDegVisual, 'linkingFuncParams', p.Results.midgetLinkingFuncParams );
    plot(regularSupportPosDegVisual,midgetFraction,'-r');
    
end

end % rgcThickness function

