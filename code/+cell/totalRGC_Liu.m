function totalRGC = totalRGC_Liu(showPlots)
% Count functions for the sum of all retinal ganglion cells
%
% Syntax:
%  totalRGC = cell.totalRGC(showPlots)
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) as a function of
%   retinal eccentricity (in degrees).
%
% Examples:
%{
    totalRGC = cell.totalRGC(true)
%}

% Handle plotting
if nargin==0
    showPlots = false;
end

% The meridians over which the calculation is to be performed
cardinalMeridianAngles = [0 180];
cardinalMeridianNames = {'nasal' 'temporal'};

%% Cell counts
% We call the rgcDisplacementMap functions, passing the measurmements taken
% from figure 3 of Liu 2017 PNAS
dataFileName = fullfile(fileparts(mfilename('fullpath')),'liuRawRGCDensity_reportedAverage.mat');

% Loop over the specified meridians
for mm = 1:length(cardinalMeridianAngles)
    
    % Set up this meridian model element
    totalRGC(mm).label = cardinalMeridianNames(mm);
    totalRGC(mm).angle = cardinalMeridianAngles(mm);
    
    % Obtain a spline fit to the cell densities
    splineFit = getSplineFitToRGCDensitySqDegVisual(cardinalMeridianAngles(mm),...
        'cardinalMeridianAngles',cardinalMeridianAngles,...
        'cardinalMeridianNames',cardinalMeridianNames,...
        'splineKnots',5,...
        'splineOrder',4,...
        'rgcDensityDataFileName',dataFileName);

    % Nan optic disc points and those points where the spline fit returns a
    % value of zero (i.e., points where we have not data regarding RGC
    % density)
    totalRGC(mm).countsDegSq =  @(posDeg) ...
        nanZeros(...
        nanOpticDiscPoints(splineFit(posDeg), posDeg, cardinalMeridianAngles(mm)) ...
        );

    if showPlots
        if mm == 1
            figure
        end
        [rgcDensitySqDegVisual, rgcNativeSupportPosDegVisual] =  ...
            loadRawRGCDensityByEccen(cardinalMeridianAngles(mm), ...
            'rgcDensityDataFileName', dataFileName);
        plot(rgcNativeSupportPosDegVisual,rgcDensitySqDegVisual,'xk')
        hold on
        plot(0:0.5:50,totalRGC(mm).countsDegSq(0:0.5:50));
    end
    
end

end

function y = nanZeros(y)
    y(y==0) = nan;
end
