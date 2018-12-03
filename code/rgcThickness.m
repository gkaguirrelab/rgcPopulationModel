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
%
% Outputs:
%
% Examples:
%{
    % Search across midget fraction linking params to obtain the best model
    % fit to tissue thickness
    myObj = @(p) rgcThickness('midgetLinkingFuncParams',p,'showPlots',false);
    x0=[12.0290    1.7850];
    [x,fval]=fmincon(myObj,x0)
%}

%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);
%p.addParameter('midgetLinkingFuncParams',[12.0290    1.7850],@isnumeric); % Dacey model
%p.addParameter('midgetLinkingFuncParams',[0.4476   48.2818],@isnumeric); % Parameters thata best fit thickness when zero tissue shrinkage is assumed
p.addParameter('midgetLinkingFuncParams',[2.1983    1.2463],@isnumeric); % Current best linking function params from the Barnett & Aguirre model
%p.addParameter('midgetLinkingFuncParams',[0.2812   27.7541],@isnumeric); % Parameters thata best fit thickness when tissue shrinkage is assumed

p.addParameter('showPlots',true,@islogical);

% parse
p.parse(varargin{:})


% Obtain the cell population components from the individual functions
amacrine = cell.amacrine(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
totalRGC = cell.totalRGC(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
midget = cell.midget(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames,p.Results.midgetLinkingFuncParams);
bistratified = cell.bistratified(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames, totalRGC);
parasol = cell.parasol(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames, totalRGC, midget, bistratified);
ipRGC = cell.ipRGC(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);

% Obtain the empirical thickness measurements
rgcLayer = layer.rgc(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
iplLayer = layer.ipl(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);

% Build a set of anonymous functions to support fitting

% Volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

% Packing density of spheres
% We model the cells as spheres. For a single sphere size, the densest
% possible packing is Keppler's limit of ~0.74. In the RGC there are
% multiple cell classes. Except in the far periphery the two most numerous
% types of cells are the midget and parasol. If smaller cells are more
% numerous than larger cells, it is possible in principle to pack small
% cells into the spaces between larger cells, and achieve a density higher
% than Kepler's limit. When the ratio of the smaller to larger cell is
% greater than ~0.41, then the smaller is too large to pack into the
% larger. This is the situation we face in this model, so we assume
% Keppler's limit for sphere packing.
spherePackDensity = 0.74048048969;


% Pick the meridian to model
meridian = 'temporal';

% Set the support to be the locations of the Curcio 2011 RGC thickness
% measurements
supportDeg = rgcLayer.supportDeg.(meridian)';

calcRGCthickness = @(xScale,yScale) ((...
    amacrine.density.fitDegSq.(meridian)(supportDeg.*xScale) .* sVol(amacrine.diameter.fitDeg(supportDeg.*xScale)) + ...
    parasol.density.fitDegSq.(meridian)(supportDeg.*xScale) .* sVol(parasol.diameter.fitDeg(supportDeg.*xScale)) + ...
    bistratified.density.fitDegSq.(meridian)(supportDeg.*xScale) .* sVol(bistratified.diameter.fitDeg(supportDeg.*xScale)) + ...
    midget.density.fitDegSq.(meridian)(supportDeg.*xScale) .* sVol(midget.diameter.fitDeg(supportDeg.*xScale)) + ...
    ipRGC.density.fitDegSq.(meridian)(supportDeg.*xScale) .* sVol(ipRGC.diameter.fitDeg(supportDeg.*xScale)) ...
    ) .* yScale) ./ spherePackDensity ./ calc_mmSqRetina_per_degSqVisual(supportDeg,180);

myObj = @(x) sqrt(sum((rgcLayer.thickMM.(meridian)' - calcRGCthickness(x(1),x(2))).^2));

[fParams, fVal] = fmincon(myObj,[1 1]);

%% Figure prep

if p.Results.showPlots
    figure
    
    % Plot counts
    subplot(2,2,1)
    %plot(supportDeg, totalRGC.density.fitDegSq.(meridian)(supportDeg))
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
    
    plot(rgcLayer.supportDeg.(meridian), rgcLayer.thickMM.(meridian), 'xk');
    hold on
    plot(supportDeg, calcRGCthickness(1,1),'-k');
    legend({'Curcio RGC measure','Adjusted model','Fixed model'});
    xlabel('eccentricity [mm retina]');
    ylabel('layer thickness [mm]]');
    
end

end % rgcThickness function

