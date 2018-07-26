function rgcThickness(varargin )
% Caclulates RGC layer thickness by reference to constituent cell classes
%
% Description:
%   We aim to create a model of the thickness of the retinal ganglion cell
%   layer and how it vaires with eccentricity, comparing our estimations to
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
%


%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('polarAngle',180,@isnumeric);
p.addParameter('cardinalMeridianAngles',[0 90 180 270],@isnumeric);
p.addParameter('cardinalMeridianNames',{'nasal' 'superior' 'temporal' 'inferior'},@iscell);

% parse
p.parse(varargin{:})


% Obtain the cell population components from the individual functions
amacrine = cell.amacrine(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
totalRGC = cell.totalRGC(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
midget = cell.midget(p.Results.cardinalMeridianAngles, p.Results.cardinalMeridianNames);
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

% Set the support to be the locations of the Curcio 2011 RGC thickness
% measurements
supportMM = rgcLayer.supportMM.temporal';

calcRGCthickness = @(xScale,yScale) ((...
    amacrine.density.fitMMSq.temporal(supportMM.*xScale) .* sVol(amacrine.diameter.fitMM(supportMM.*xScale)) + ...
    parasol.density.fitMMSq.temporal(supportMM.*xScale) .* sVol(parasol.diameter.fitMM(supportMM.*xScale)) + ...
    bistratified.density.fitMMSq.temporal(supportMM.*xScale) .* sVol(bistratified.diameter.fitMM(supportMM.*xScale)) + ...
    midget.density.fitMMSq.temporal(supportMM.*xScale) .* sVol(midget.diameter.fitMM(supportMM.*xScale)) + ...
    ipRGC.density.fitMMSq.temporal(supportMM.*xScale) .* sVol(ipRGC.diameter.fitMM(supportMM.*xScale)) ...
    ) .* yScale) ./ spherePackDensity;

myObj = @(x) sqrt(sum((rgcLayer.thickMM.temporal' - calcRGCthickness(x(1),x(2))).^2));

[fParams, fVal] = fmincon(myObj,[1 1]);

figure
plot(supportMM,rgcLayer.thickMM.temporal,'-k');
hold on
plot(supportMM,calcRGCthickness(fParams(1),fParams(2)),'*r');


%% Figure prep
figure

% Plot counts
subplot(1,3,1)
plot(supportMM, totalRGC.density.fitMMSq.temporal(supportMM))
hold on
plot(supportMM, midget.density.fitMMSq.temporal(supportMM))
plot(supportMM, parasol.density.fitMMSq.temporal(supportMM))
plot(supportMM, bistratified.density.fitMMSq.temporal(supportMM))
plot(supportMM, amacrine.density.fitMMSq.temporal(supportMM))
plot(supportMM, ipRGC.density.fitMMSq.temporal(supportMM))
plot(supportMM, midget.density.fitMMSq.temporal(supportMM) + parasol.density.fitMMSq.temporal(supportMM) + bistratified.density.fitMMSq.temporal(supportMM) +ipRGC.density.fitMMSq.temporal(supportMM) + amacrine.density.fitMMSq.temporal(supportMM),'xr')
xlabel('eccentricity [mm retina]');
ylabel('density [counts / sq mm]');
legend({'Curcio totalRGC','midget','parasol','bistratified','amacrine','ipRGC','model total all cells'});

% Plot cell volumes


subplot(1,3,2)
plot(supportMM, sVol(midget.diameter.fitMM(supportMM)))
hold on
plot(supportMM, sVol(parasol.diameter.fitMM(supportMM)))
plot(supportMM, sVol(bistratified.diameter.fitMM(supportMM)))
plot(supportMM, sVol(amacrine.diameter.fitMM(supportMM)))
plot(supportMM, sVol(ipRGC.diameter.fitMM(supportMM)))
xlabel('eccentricity [mm retina]');
ylabel('individual cell volume [mm^3]');
legend({'midget','parasol','bistratified','amacrine','ipRGC'});


% Plot thickness
subplot(1,3,3)

plot(supportMM, calcRGCthickness(fParams(1),fParams(2)));
hold on
plot(rgcLayer.supportMM.temporal, rgcLayer.thickMM.temporal, '*r');
legend({'Model thickness','Curcio RGC measure'});
xlabel('eccentricity [mm retina]');
ylabel('layer thickness [mm]]');




end % rgcThickness function


%% LOCAL FUNCTIONS
