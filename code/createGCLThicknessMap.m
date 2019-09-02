function gclThicknessMap = createGCLThicknessMap(varargin)
% Caclulates RGC layer thickness by reference to constituent cell classes
%
% Description:

%
% Inputs:
%
% Optional key / value pairs:
%  'packingDensity'       - Scalar. We model the cells as spheres. For a

% Outputs:
%
% Examples:
%{
    gclThicknessMap = createGCLThicknessMap();
    imagesc(gclThicknessMap);
    axis square
%}

%% input parser
p = inputParser;

% Optional analysis params
p.addParameter('octDataFileName', ...
    fullfile(tbLocateProject('rgcPopulationModel','verbose',false),'data','rgcIplThicknessMap.mat'),...
    @ischar);

% parse
p.parse(varargin{:})

% Load the GCL+IPL map
load(p.Results.octDataFileName,'rgcIplThicknessMap');
octRadialDegreesVisualExtent = 15;

% Convert from microns to mm
rgcIplThicknessMap = rgcIplThicknessMap./1000;

scaleDown = 12;
supportDeg = linspace(-octRadialDegreesVisualExtent,octRadialDegreesVisualExtent,size(rgcIplThicknessMap,1)/scaleDown);
downSampMap = imresize(rgcIplThicknessMap,1/scaleDown);
[xo,yo,zo]=prepareSurfaceData(supportDeg,supportDeg,downSampMap);
sf = fit([xo, yo],zo,'thinplateinterp');
fitHandle = @(x,y) sf(x,y);
yMax = @(x) sqrt(octRadialDegreesVisualExtent^2 - x.^2);
yMin = @(x) -sqrt(octRadialDegreesVisualExtent^2 - x.^2);
totalRGCIPLthickness = integral2(fitHandle,-octRadialDegreesVisualExtent,octRadialDegreesVisualExtent,yMin,yMax);

surf(supportDeg,supportDeg,downSampMap)
hold on
[xi,yi]=meshgrid(-octRadialDegreesVisualExtent:0.25:octRadialDegreesVisualExtent,-octRadialDegreesVisualExtent:0.25:octRadialDegreesVisualExtent);
inRangeIdx = sqrt(xi.^2+yi.^2)<=octRadialDegreesVisualExtent;
plot3(xi(inRangeIdx),yi(inRangeIdx),sf(xi(inRangeIdx),yi(inRangeIdx)),'.r')

% Create a version of the original GCL+IPL map which has interpolated
% values for the missing points
supportDeg = linspace(-octRadialDegreesVisualExtent,octRadialDegreesVisualExtent,size(rgcIplThicknessMap,1));
[xi,yi]=meshgrid(supportDeg,supportDeg);
inRangeIdx = sqrt(xi.^2+yi.^2)<=octRadialDegreesVisualExtent;
nanIdx = isnan(rgcIplThicknessMap);
fillInPoints = logical(inRangeIdx.*nanIdx);

rgcIplThicknessMapFixed = nan(size(rgcIplThicknessMap));
rgcIplThicknessMapFixed(inRangeIdx)=rgcIplThicknessMap(inRangeIdx);
rgcIplThicknessMapFixed(fillInPoints) = sf(xi(fillInPoints),yi(fillInPoints));

% Convert from an image to a polar representation
polarMap = convertImageMapToPolarMap(rgcIplThicknessMapFixed);
gclThicknessMapPolar = nan(size(polarMap));

ratioFuncByThickness = rgcLayerProportion('referenceEccenDegVisual',octRadialDegreesVisualExtent);

% Obtain the cumulative thickness along each row
profileSupport = linspace(0,octRadialDegreesVisualExtent,size(polarMap,2));
for ii=1:size(polarMap,1)
    gclIPL_profile = polarMap(ii,:);
    gcl_profile = gclIPL_profile.*ratioFuncByThickness(profileSupport,gclIPL_profile,totalRGCIPLthickness);
    gclThicknessMapPolar(ii,:) = gcl_profile;
end


gclThicknessMap = convertPolarMapToImageMap(gclThicknessMapPolar);


end

