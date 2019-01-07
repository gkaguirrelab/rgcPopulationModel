function fVal=rgcThickness(varargin )
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
    % Search across midget fraction linking params to obtain the best model
    % fit to tissue thickness
    myObj = @(p) rgcThickness('midgetLinkingFuncParams',p,'showPlots',false,'forceRecalculate',false);
    x0=[12.0290    1.7850];
    ub=[22 3];
    lb=[6 1];
    [x,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
%}
%{
    % Search across packing density to obtain the best model % fit to
    tissue thickness
    myObj = @(p) rgcThickness('packingDensity',p,'showPlots',false,'forceRecalculate',false);
    x0=0.74;
    [x,fval]=fmincon(myObj,x0,[],[])
%}
%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);
%p.addParameter('midgetLinkingFuncParams',[12.0290    1.7850],@isnumeric); % Dacey params
p.addParameter('midgetLinkingFuncParams',[6.0007    2.0455],@isnumeric); % Best fit to the OCT data
%p.addParameter('midgetLinkingFuncParams',[2.2, 1.25],@isnumeric); % Barnett Aguirre values
%p.addParameter('packingDensity',0.74048048969,@isscalar);
p.addParameter('packingDensity',0.5428,@isscalar);  % Best fit to the OCT data

p.addParameter('forceRecalculate',false,@islogical);


p.addParameter('showPlots',true,@islogical);
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
persistent rgcOCTMm
if isempty(rgcOCTMm) || p.Results.forceRecalculate
    % Load the OCT data
    dataLoad = load(p.Results.octDataFileName);
    rgcIplThicknessMap = dataLoad.rgcIplThicknessMap;
    octRadialDegreesVisualExtent = 15;
    clear dataLoad
    
    % Extract the horizontal meridian
    rgciplOCTthickness = rgcIplThicknessMap(round(size(rgcIplThicknessMap,1)/2),:);
    
    % Construct a data structure to hold the OCT thickness values. Start with
    % the support in visual degrees
    rgcOCTMm.supportDeg.temporal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*octRadialDegreesVisualExtent;
    rgcOCTMm.supportDeg.nasal = ((1:round(length(rgciplOCTthickness)/2))./round(length(rgciplOCTthickness)/2)).*octRadialDegreesVisualExtent;
    
    % Now obtain the ratio of RGC thickness to RGC+IPL thickness
    thicknessRatioTemporal = rgcIplThicknessRatio( rgcOCTMm.supportDeg.temporal, 'forceRecalculate', true)';
    thicknessRatioNasal = rgcIplThicknessRatio( rgcOCTMm.supportDeg.nasal, 'forceRecalculate', true )';
    
    % Calculate and store the thickness of the RGC layer in mm
    rgcOCTMm.thickMM.temporal = fliplr(rgciplOCTthickness(1:round(length(rgciplOCTthickness)/2))).*thicknessRatioTemporal./1000;
    rgcOCTMm.thickMM.nasal = rgciplOCTthickness(round(length(rgciplOCTthickness)/2)+1:end).*thicknessRatioNasal./1000;
    
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

% Set the support to be the locations of the Curcio 2011 RGC thickness
% measurements
supportDeg = (0:1:10)';

for mm=1:2
    calcRGCthickness.(cardinalMeridianNames{mm}) = ((...
        amacrine.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(amacrine.diameter.fitDeg(supportDeg)) + ...
        parasol.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(parasol.diameter.fitDeg(supportDeg)) + ...
        bistratified.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(bistratified.diameter.fitDeg(supportDeg)) + ...
        midget.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(midget.diameter.fitDeg(supportDeg)) + ...
        ipRGC.density.fitDegSq.(cardinalMeridianNames{mm})(supportDeg) .* sVol(ipRGC.diameter.fitDeg(supportDeg)) ...
        )) ./ p.Results.packingDensity ./ calc_mmSqRetina_per_degSqVisual(supportDeg,180);
end

fVal = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDeg') - calcRGCthickness.nasal).^2)) + ...
    sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDeg') - calcRGCthickness.temporal).^2));

fVal = sqrt(sum((rgcOCTMm.fitDeg.nasal(supportDeg')./max(rgcOCTMm.fitDeg.nasal(supportDeg')) - calcRGCthickness.nasal./max(calcRGCthickness.nasal)).^2)) + ...
    sqrt(sum((rgcOCTMm.fitDeg.temporal(supportDeg')./max(rgcOCTMm.fitDeg.temporal(supportDeg')) - calcRGCthickness.temporal./max(calcRGCthickness.temporal)).^2));

fVal

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
end

end % rgcThickness function

