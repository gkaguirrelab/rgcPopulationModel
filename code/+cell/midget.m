function midget = midget( cardinalMeridianAngles, cardinalMeridianNames )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%% Midget RGCs

% Obtain a spline fit to RGC density values from Curcio & Allen
for mm = 1:length(cardinalMeridianAngles)
    midget.density.fitMMSq.(cardinalMeridianNames{mm}) = createMidgetDensityFunc(cardinalMeridianAngles(mm));
end

% Midget cell body sizes Liu and colleagues 2017 report midget soma sizes
% as a function of eccentricity as measured by AO-OCT:
%   Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other
%   transparent neurons in the living human retina." Proceedings of the
%   National Academy of Sciences (2017): 201711734.
%
% Dacey 1993 (Table 1) reports a larget size of 17.9 microns.
%
%   Dacey, Dennis M. "Morphology of a small-field bistratified ganglion
%   cell type in the macaque and human retina." Visual neuroscience 10.6
%   (1993): 1081-1098.

% Support in the source data is in degrees of visual field along the
% temporal retina
midget.diameter.supportDeg = [(1.5+3)/2, (3+4.5)/2, (6+7.5)/2, (8+9.5)/2, (12+13.5)/2];

% We convert from visual degrees back to retinal mm.
midget.diameter.supportMM = convert_degVisual_to_mmRetina(midget.diameter.supportDeg, 180);
midget.diameter.sizeMM = [0.0115, 0.0113, 0.0114, 0.0118, 0.01315];

% Obtain a spline fit
midget.diameter.fitMM = fit(midget.diameter.supportMM', midget.diameter.sizeMM', 'poly1');

end



function midgetDensityFitMMSq = createMidgetDensityFunc(meridianAngle)

% The Barnett & Aguirre 2018 Linking Function parameters for each of the
% cardinal meridians (0, 90, 180, 270).
midgetFracLinkingParams = [...
   4.3286    1.6160; ...
   4.2433    1.5842; ...
   4.2433    1.5842; ...
   4.2433    1.6160];

% Define a support vector in retinal degrees
regularSupportPosDegRetina = 0:0.01:30;

% Obtain the spline fit function to total RGC density
fitRGCDensitySqDegRetina = getSplineFitToRGCDensitySqDegRetina(meridianAngle);

% Define a variable with RGC density over regular support and zero values
% at the optic disc positions
RGCDensityOverRegularSupport = ...
    zeroOpticDiscPoints(fitRGCDensitySqDegRetina(regularSupportPosDegRetina),regularSupportPosDegRetina, meridianAngle);

% Define the mRGC density function using the Barnett & Aguirre 2018 linking
% function parameters
mRGCDensityOverRegularSupport = transformRGCToMidgetRGCDensityDacey(regularSupportPosDegRetina,RGCDensityOverRegularSupport,...
    'linkingFuncParams',midgetFracLinkingParams(meridianAngle/90+1,:));

% Perform a spline fit to the mRGC density expressed in mm and mm sq units
midgetDensityFitMMSq = fit(convert_degRetina_to_mmRetina(regularSupportPosDegRetina)', (mRGCDensityOverRegularSupport.*calc_degSqRetina_per_mmSqRetina())', 'smoothingspline');

end

function vectorOut = zeroOpticDiscPoints(vectorIn, regularSupportPosDegRetina, polarAngle)
opticDiscIndices = findOpticDiscPositions(regularSupportPosDegRetina, polarAngle);
vectorOut = vectorIn;
vectorOut(opticDiscIndices) = 0;
end


