function fVal=modelGCLayerThickness(varargin )
% Caclulates GC layer thickness by reference to constituent cell classes
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
    modelGCLayerThickness()
%}
%{
    % Search across model params to fit the data
    myObj = @(p) modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',false,'forceRecalculate',false);
    x0=[14.9920    0.6998    0.5990    0.9000  0.6122];
    ub=[15 0.7 0.6 1.0 0.7];
    lb=[1 0.1 0.2 0.9 0.5];
    [p,fval]=fmincon(myObj,x0,[],[],[],[],lb,ub)
    modelGCLayerThickness('midgetLinkingFuncParams',p(1:4),'packingDensity',p(5),'showPlots',true,'forceRecalculate',false,'objectiveType','all');
%}

%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('midgetLinkingFuncParams',[14.9920    0.6998    0.5990    0.9000],@isnumeric); % Best fit to the OCT data
p.addParameter('packingDensity',0.6199,@isscalar);  % Best fit to the OCT data
p.addParameter('forceRecalculate',true,@islogical);
p.addParameter('showPlots',true,@islogical);
p.addParameter('outDir','~/Desktop/KaraCloud_VSS2019_figs',@ischar);

% parse
p.parse(varargin{:})

% Obtain the ganglion cell thickness data
data = prepareGCThickData();

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
midget = cell.midget(totalRGC, p.Results.midgetLinkingFuncParams);
parasol = cell.parasol(totalRGC, midget, bistratified);

% An anonymous function for the volume of a sphere given diameter
sVol = @(d) 4/3*pi*(d./2).^3;

% Loop over the meridians provided in the data. Build the model and obtain
% the error.
fVal = 0;
for mm=1:length(data)
    
    % Prepare this entry of the model
    model(mm).label = data(mm).label;
    model(mm).angle = data(mm).angle;
    
    % Find the matching cell component entry
    idx = strcmp([totalRGC(:).label],data(mm).label);

    % The model thickness is the sum of the layer components
    model(mm).thickMM = @(x) ((...
        amacrine(idx).countsDegSq(x) .* sVol(amacrine(idx).diameter(x)) + ...
        parasol(idx).countsDegSq(x) .* sVol(parasol(idx).diameter(x)) + ...
        bistratified(idx).countsDegSq(x) .* sVol(bistratified(idx).diameter(x)) + ...
        midget(idx).countsDegSq(x) .* sVol(midget(idx).diameter(x)) + ...
        ipRGC(idx).countsDegSq(x) .* sVol(ipRGC(idx).diameter(x)) ...
        )) ./ p.Results.packingDensity ./ calc_mmSqRetina_per_degSqVisual(x',totalRGC(mm).angle);

    % Evaluate the model at the data locations
    dataModelDif = data(mm).thickMM - model(mm).thickMM(data(mm).supportDeg);
    fVal = fVal + sqrt(nansum(dataModelDif.^2));
    
end


end 

