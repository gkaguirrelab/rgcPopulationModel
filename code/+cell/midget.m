function midget = midget( totalRGC, midgetLinkingFuncParams, showPlots )
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

% Handle plotting and missing midgetLinkingFuncParams
if nargin==1
    midgetLinkingFuncParams = [12.0290, 0.4320];
    showPlots = false;
end

if nargin==2
    showPlots = false;
end

if showPlots
    figure
end

% Define a maximum eccentricity of the model
maxEccenDeg = 50;


%% Cell counts
% The midgetLinkingFuncParams are used to construct a function that
% converts the values within totalRGC into the midget fraction

% Define a support vector in visual degrees
supportDeg = 0:0.01:maxEccenDeg;

% Loop over the specified meridians
for mm = 1:length(totalRGC)
    
    % Obtain the total RGC density for this support
    countsDegSqTotal = totalRGC(mm).countsDegSq(supportDeg);
    
    % Set the optic disc nans to zero
    countsDegSqTotal(isnan(countsDegSqTotal)) = 0;
    
    % Define a function that returns the midget density, conditioned on the
    % property of the parameters
    switch length(midgetLinkingFuncParams)
        case 1
            countsDegSq = transformRGCToMidgetRGCDensityDrasdo(supportDeg,countsDegSqTotal);
        case 2
            countsDegSq = transformRGCToMidgetRGCDensityDacey(supportDeg,countsDegSqTotal,...
                'linkingFuncParams',midgetLinkingFuncParams);
        case 3
            countsDegSq = transformRGCToMidgetRGCDensityDacey(supportDeg,countsDegSqTotal,...
                'linkingFuncParams',midgetLinkingFuncParams(1:2),'minMidgetFractionRatio',midgetLinkingFuncParams(3));
        case 4
            countsDegSq = transformRGCToMidgetRGCDensityDacey(supportDeg,countsDegSqTotal,...
                'linkingFuncParams',midgetLinkingFuncParams(1:2),'minMidgetFractionRatio',midgetLinkingFuncParams(3),'maxMidgetFractionRatio',midgetLinkingFuncParams(4));
    end
    
    % Obtain a fit to the cell densities
    smoothFit = fit(supportDeg', countsDegSq', 'cubicinterp');

    % Set up this meridian model element
    midget(mm).label = totalRGC(mm).label;
    midget(mm).angle = totalRGC(mm).angle;
    
    % Nan optic disc points and save the anonymous function
    midget(mm).countsDegSq = @(posDeg) ...        
        nanOpticDiscPoints(smoothFit(posDeg), posDeg, totalRGC(mm).angle);

    % Plot the fit
    if showPlots
        plot(supportDeg,midget(mm).countsDegSq(supportDeg));
        hold on
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

if showPlots
    figure
end

% Loop over the specified meridians
for mm = 1:length(totalRGC)

    % Support in the source data is in degrees of visual field along the
    % temporal retina
    supportDeg = [(1.5+3)/2, (3+4.5)/2, (6+7.5)/2, (8+9.5)/2, (12+13.5)/2];
    sizeMM = [0.0115, 0.0113, 0.0114, 0.0118, 0.01315];

    % Obtain the fit and save. We find that an exponential does best with
    % the impoverished set of measurements we have
    fx = @(a,b,c,x) (a.*x).^b+c;
    midget(mm).diameter = fit(supportDeg', sizeMM',...
        fx,'StartPoint', [0.001 1 0.0115], ...
        'Lower', [0 1 0]);

    % Plot the fit
    if showPlots
        plot(0:0.01:maxEccenDeg,midget(mm).diameter(0:0.01:maxEccenDeg));
        hold on
        plot(supportDeg,sizeMM,'*');
    end

end

end


