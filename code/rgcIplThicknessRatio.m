function [ rgcProportionThickness ] = rgcIplThicknessRatio( regularSupportPosDegVisual, varargin )
% Returns the proportion of RGC+IPL layer thickness on OCT that is RGC
%
% Description:
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
%   the basis of Figure 7.
%
%   These measurements were made on eyes from older (40-92) people, and
%   there are various changes to tissue thickness that can happen as a
%   consequence of fixation. While the absolute thickness values may be
%   inaccurate, we are willing to assume that the ratio of RGC and IPL
%   layer thickness is accurate.
%
% Inputs:
%   regularSupportPosDegVisual - A 1 x p vector that contains the
%                           eccentricity in visual degrees at which the
%                           model was evaluated along each meridian
%
% Optional key/value pairs:
%  'forceRecalculate'     - When set to true forces the routine to
%                           re-calculate the thickness ratio function from
%                           the source data.
%
% Outputs:
%   rgcProportionThickness - a 1xp vector that provides the proportion of
%                           thickness at each eccentricity location which
%                           is RGC out of the RGC+IPL layer.
%
% Examples:
%{
    supportDeg = 0:0.1:30;
    rgcProportionThickness = rgcIplThicknessRatio(supportDeg);
    plot(supportDeg,rgcProportionThickness)
%}


%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDegVisual',@isnumeric);

% optional key/value pairs
p.addParameter('forceRecalculate',true,@islogical);


% parse
p.parse(regularSupportPosDegVisual, varargin{:})

% Define the persistent variable ratioFit. This way, we do not need to
% re-compute this spline fit result each time we encounter this file
persistent ratioFit

if isempty(ratioFit) || p.Results.forceRecalculate
    
    % Obtain the Curcio 2011 empirical measurements of RGC and IPL layer
    % thickness from histology
    rgcLayer = layer.rgc();
    iplLayer = layer.ipl();
    
    % Across regular support, obtain the proportion of rgc layer thickness and
    % then fit this with a cubic spline. This manuever is needed to prevent
    % strange value occurring at very small eccentricities
    localSupport = 0.5:0.1:15;
    rgcProportionThickness = mean([...
        (rgcLayer.fitDeg.temporal(localSupport)./(rgcLayer.fitDeg.temporal(localSupport)+iplLayer.fitDeg.temporal(localSupport)))'; ...
        (rgcLayer.fitDeg.nasal(localSupport)./(rgcLayer.fitDeg.nasal(localSupport)+iplLayer.fitDeg.nasal(localSupport)))']);
    ratioFit = fit([0 localSupport]', [0 rgcProportionThickness]', 'cubicspline');
    
end

% Now evaluate the ratio function at the passed support positions
rgcProportionThickness = ratioFit(regularSupportPosDegVisual);

end % function


