function ratioFuncByThickness = rgcLayerProportion( varargin )
% A model for the proportion of RGC+IPL layer thickness on OCT that is RGC
%
% Syntax:
%  ratioFunc = rgcLayerProportion
%
% Description:
%   This routine creates a function that models the relationship between
%   the thickness of the RGC+IPL layer of the retina and the proportion of
%   that thickness that is the RGC component.
%
%   Curcio and colleagues measured on histological sections the thickness
%   of the RGC and IPL layers of the retina in a set of ex-vivo eyes. These
%   measurements were made at positions along the nasal and temporal
%   meridian, expressed in mm from the fovea:
%
%       Curcio, Christine A., et al. "Human chorioretinal layer thicknesses
%       measured in macula-wide, high-resolution histologic sections."
%       Investigative ophthalmology & visual science 52.7 (2011):
%       3943-3954.
%
%   Christine provided us with a spreadsheet with the measurements that are
%   the basis of Figure 7 in her paper. These values are stored in the
%   routines layer.ipl and layer.rgc.
%
%   These measurements were made on eyes from older (40-92) people, and
%   there are various changes to tissue thickness that can happen as a
%   consequence of fixation. While the absolute thickness values may be
%   inaccurate, we are willing to assume that the ratio of RGC and IPL
%   layer thickness is accurate.
%
%   The model 
%
% Inputs:
%   none
%
% Optional key/value pairs:b
%  'referenceEccenDegVisual' - Scalar. The reference eccentricity for the
%                           proportion of the cumulative layer thickness.
%                           The proportion function will have a value of
%                           unity at this point.
%  'showPlots'            - Logical.
%
% Outputs:
%   ratioFuncByThickness  - A function handle that takes as input the
%                           1xn vectors regularSupportVisDeg and thickness
%                           and returns the 1xn vector that is the
%                           proportion of the thickness at each support
%                           position that is RGC.
%
% Examples:
%{
    rgcLayerProportion('showPlots',true);
%}


%% Parse input and define variables
p = inputParser;

% optional key/value pairs
p.addParameter('referenceEccenDegVisual',15,@isscalar);
p.addParameter('showPlots',false,@islogical);

% parse
p.parse(varargin{:})


% Obtain the Curcio 2011 empirical measurements of RGC and IPL layer
% thickness from histology
rgcLayer = layer.rgc();
iplLayer = layer.ipl();

% Regular support for position om the retina
regularSupportVisDeg = 0.5:0.01:15;

% Regular support for thickness cumulative sum
regularSupportCumSum = 0.0005:0.0005:1.5;

% The meridians for which we have data (should be temporal and nasal)
meridians = fieldnames(rgcLayer.fitDeg);

% Create a variable to hold the obtained relationship across meridians
ratioRelationship = zeros(size(regularSupportCumSum));

% Plot
if p.Results.showPlots
    figure
end

% Loop over the meridians to obtain the cumulative sum thickness at the
% reference eccentricity
[~, referenceIdx] = min(abs(regularSupportVisDeg-p.Results.referenceEccenDegVisual));
for mm = 1:length(meridians)
    
    % The RGC+IPL layer thickness across regular support
    thickRGCIPL = (rgcLayer.fitDeg.(meridians{mm})(regularSupportVisDeg) + ...
        iplLayer.fitDeg.(meridians{mm})(regularSupportVisDeg))';
    thicknessCumSum = calcRingCumulative(regularSupportVisDeg,thickRGCIPL);
    
    referenceCumSumThick(mm) = thicknessCumSum(referenceIdx);
end
referenceCumSumThick = mean(referenceCumSumThick);

% Now loop again and obtain the cumulative sum thickness relative to the
% reference.
for mm = 1:length(meridians)
    
    % Obtain the cumulative sum of total thickness of the RGC+IPL layers
    % across regular support, relative to the reference eccentricity
    thickRGCIPL = (rgcLayer.fitDeg.(meridians{mm})(regularSupportVisDeg) + ...
        iplLayer.fitDeg.(meridians{mm})(regularSupportVisDeg))';
    thicknessCumSum = ...
        relativeCumulativeThickness(regularSupportVisDeg, thickRGCIPL, referenceCumSumThick);
    
    % Obtain the ratio [ RGC / (RGC+IPL) ] across regular support
    rgcThickRatio = rgcLayer.fitDeg.(meridians{mm})(regularSupportVisDeg)' ./ thickRGCIPL;
    
    % Plot
    if p.Results.showPlots
        subplot(1,4,1)
        plot(regularSupportVisDeg,rgcThickRatio)
        xlabel('Eccentricity visual degrees');
        ylabel('Proportion of RGC+IPL layer that is RGC');
        title('RGC proportion of RGCIPL');
        hold on
        subplot(1,4,2)
        plot(regularSupportVisDeg,thicknessCumSum)
        xlabel('Eccentricity visual degrees');
        ylabel('Relative cumulative sum thickness');
        title('cumulative sum thickness RGCIPL');
        hold on
    end
    
    % Obtain a spline fit to the relationship between relative
    % cumulative thickness and the thickness ratio
    ratioRelationshipByMeridianFit = spline(thicknessCumSum,rgcThickRatio);
    
    % Add this fitted relationship to the variable that accumulates across
    % loops over meridians
    ratioRelationship = ratioRelationship + ...
        ppval(ratioRelationshipByMeridianFit,regularSupportCumSum);
    
    % Plot
    if p.Results.showPlots
        subplot(1,4,3)
        plot(thicknessCumSum,rgcThickRatio)
        hold on
    end
end

% Obtain the average relationship across meridians
ratioRelationship = ratioRelationship ./ length(meridians);

% Create a function handle that provides the ratio relationship for any
% particular cumulative sum thickness value
ratioFuncByCumThickness = @(x) ppval(spline(regularSupportCumSum,ratioRelationship),x);

% And now a function to return that provides the thickness of the RGC
% portion of the RGC+IPL layer given that total thickness over regular
% support.
ratioFuncByThickness = @(regularSupportVisDeg, thickness) ratioFuncByCumThickness(relativeCumulativeThickness(regularSupportVisDeg, thickness, referenceCumSumThick));

% Plot
if p.Results.showPlots
    subplot(1,4,3)
    plot(regularSupportCumSum,ratioRelationship,'-k')
    legend([meridians;'average']);
    xlabel('Relative cumulative thickness of RGC and IPL layers');
    ylabel('Proportion of RGC+IPL layer that is RGC');
    title('RGC proportion of RGCIPL');
    
    subplot(1,4,4)
    hold on
    expandedRegularSupportVisDeg = 0:0.1:20;
    for mm=1:length(meridians)
        thickRGCIPL = (rgcLayer.fitDeg.(meridians{mm})(expandedRegularSupportVisDeg) + ...
            iplLayer.fitDeg.(meridians{mm})(expandedRegularSupportVisDeg))';
        plot(expandedRegularSupportVisDeg,ratioFuncByThickness(expandedRegularSupportVisDeg,thickRGCIPL))
    end
    xlabel('Eccentricity visual degrees');
    ylabel('Proportion of RGC+IPL layer that is RGC');
    title('RGC proportion of RGCIPL');
    legend(meridians);
end


end % function


%% LOCAL

function thicknessCumSum = relativeCumulativeThickness(regularSupportVisDeg, thickness, referenceThickness)

% Obtain the cumulative sum of the thickness across regular support
thicknessCumSum = calcRingCumulative(regularSupportVisDeg,thickness);

% Express the cumulative thickness relative to the reference thickness
thicknessCumSum = thicknessCumSum ./ referenceThickness;

end