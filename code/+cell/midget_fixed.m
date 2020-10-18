function midget = midget_fixed( totalRGC, midgetLinkingFuncParams, cellSizeParams, showPlots )
% Size and count functions for the midget RGC class
%
% Syntax:
%  midget = cell.midget( totalRGC, midgetLinkingFuncParams, showPlots )
%
% Description:
%   Returns an array of structures, where each structure has a handle to a
%   function that returns cell counts (per degree squared) and cell
%   diameter (in mm) as a function of retinal eccentricity (in degrees).
%
% Examples:
%{
	totalRGC = cell.totalRGC();
	midget = cell.midget( totalRGC, [12, 0.7, 0.41, 0.95], true );
%}

% Handle plotting and missing midgetLinkingFuncParams
if nargin==1
    midgetLinkingFuncParams = [6.2665,21.9499,0.4500,0.9500];
    showPlots = false;
end

if nargin==3
    showPlots = false;
end

% Define a maximum eccentricity and support of the model
maxEccenDeg = 50;
supportDeg = 0:0.01:maxEccenDeg;


% Obtain the midget fraction logistic function
[~, logisticFunc] = daceyMidgetFractionByEccenDegVisual();
midgetFraction = logisticFunc(...
    midgetLinkingFuncParams(1),...
    midgetLinkingFuncParams(2),...
    midgetLinkingFuncParams(3),...
    midgetLinkingFuncParams(4),...
    supportDeg);


%% Cell counts
% The midgetLinkingFuncParams are used to construct a function that
% converts the values within totalRGC into the midget fraction

% Define a support vector in visual degrees

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    % Obtain the total RGC density for this support
    countsDegSqTotal = totalRGC(mm).countsDegSq(supportDeg);
    
    % Set the optic disc nans to zero
    opticDiscPoints = isnan(countsDegSqTotal);
    
    % Obtain the midget counts
    countsDegSq = countsDegSqTotal .* midgetFraction;
    
    % Obtain a fit to the cell densities
    smoothFit = fit(supportDeg(~opticDiscPoints)', countsDegSq(~opticDiscPoints)', 'cubicinterp');
    
    % Set up this meridian model element
    midget(mm).label = totalRGC(mm).label;
    midget(mm).angle = totalRGC(mm).angle;
    
    % Nan optic disc points and save the anonymous function
    midget(mm).countsDegSq = @(posDeg) ...
        nanOpticDiscPoints(smoothFit(posDeg), posDeg, totalRGC(mm).angle);
    
    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        subplot(ceil(length(totalRGC)/2),ceil(length(totalRGC)/2),mm)
        plot(supportDeg,countsDegSqTotal,'-k');
        hold on
        plot(supportDeg,midget(mm).countsDegSq(supportDeg),'-r');
        title(totalRGC(mm).label);
    end
    
end


%% Cell diameters
% Midget cell body sizes Liu and colleagues 2017 report midget soma sizes
% as a function of eccentricity as measured by AO-OCT:
%
%   Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other
%   transparent neurons in the living human retina." Proceedings of the
%   National Academy of Sciences (2017): 201711734.
%
% Conflicting measurements are found in:
%
%   Dacey, Dennis M. "The mosaic of midget ganglion cells in the human
%   retina." Journal of Neuroscience 13.12 (1993): 5334-5355.
%
% Where the midget soma size at 1.5 mm eccentricity is reported as 17.4 or
% 18.6 microns (for the on and off midgets).

% Support in the source data is in degrees of visual field along the
% temporal retina
supportDeg = [(1.5+3)/2, (3+4.5)/2, (6+7.5)/2, (8+9.5)/2, (12+13.5)/2];
sizeMM = [0.0113, 0.0113, 0.0114, 0.0118, 0.01315];
meanSize = mean(sizeMM);
meanSupport = mean(supportDeg);

% Model the size as mean with proportional growth slope
myCellSize = @(x) (meanSize + meanSize.*((x-meanSupport).*cellSizeParams(1)+cellSizeParams(2)) )';

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    midget(mm).diameter = myCellSize;
    
    % Plot the fit
    if showPlots
        if mm == 1
            figure
        end
        plot(supportDeg,sizeMM,'x');
        hold on
        plot(0:0.01:maxEccenDeg,midget(mm).diameter(0:0.01:maxEccenDeg));
        title('modeled midget soma size')
    end
    
end

end


