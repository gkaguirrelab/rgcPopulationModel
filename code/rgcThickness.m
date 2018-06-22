function rgcThickness(varargin )
% Caclulates RGC layer thickness by reference to constituent cell classes
%
% Description:
%   Here goes some text.
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


%% Total RGC density function
% Create a function to return total RGC densities per square mm of retina
% at a given eccentricity position in mm retina.
% To do so, we call the rgcDisplacementMap functions which have a
% representation of the Curcio & Allen 1990 RGC density results. The
% rgcDisplacementMap toolbox is referenced in eccentricity units of degrees
% retina, and provides densities in square degrees retina. We convert those
% values here to mm and square mm.
for mm = 1:length(p.Results.cardinalMeridianAngles)
    tmpFit = getSplineFitToRGCDensitySqDegRetina(p.Results.cardinalMeridianAngles(mm));
    totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = @(posMMretina) tmpFit(convert_mmRetina_to_degRetina(posMMretina))'.*calc_degSqRetina_per_mmSqRetina();
end


%% Displaced amacrine cells
% Curcio & Allen 1990 distinguished displaced amacrine cells from retinal ganglion cells through imaging
% and evalutating their morphology, and determined their soma size and densities at eccentricities
% across the human retina. Amacrine cell densities are averages across four meridians.
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.

for mm = 1:length(p.Results.cardinalMeridianAngles)
    % Data from Curcio & Allen 1990, Figure 10, used for all meridians
    amacrine.density.supportMM.(p.Results.cardinalMeridianNames{mm}) = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
    amacrine.density.countsMMSq.(p.Results.cardinalMeridianNames{mm}) = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];
    % Obtain a spline fit to the amacrine densities
    amacrine.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit(amacrine.density.supportMM.(p.Results.cardinalMeridianNames{mm})', amacrine.density.countsMMSq.(p.Results.cardinalMeridianNames{mm})', 'smoothingspline');
end

% Amacrine cell body sizes Data from Curcio & Allen 1990, Figure 3 The
% support for this is the frequency it is measured in, I assumed
% percentages but it adds up to well over 100%, ask about x & y axes
amacrine.diameter.supportMM = [1, 3.95, 11.95, 55, 34.1, 11, 13, 6, 4.95, 0.8, 0.9];
amacrine.diameter.sizeMM = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17];
amacrine.diameter.fitMM = @(x) 0.01;

%% Parasol RGCs
% Average parasol cell densities across four meridians measured in
% cells/square mm from six macaque retinas
% 
%   Silveira, L. C. L., and V. H. Perry. "The topography of magnocellular
%   projecting ganglion cells (M-ganglion cells) in the primate retina."
%   Neuroscience 40.1 (1991): 217-237.

% % Data from Silviera et al. 1991, Figures 16
% parasol.density.supportMM.temporal = 1:1:22;
% parasol.density.countsMMSq.temporal = [2150, 1550, 900, 460, 300, 220, 170, 110, 70, 56, 40, 20, 10, nan, nan, nan, nan, nan, nan, nan, nan, nan];
% parasol.density.supportMM.nasal = 1:1:22;
% parasol.density.countsMMSq.nasal = [3350, 2560, nan, nan, nan, 520, 485, 415, 400, 380, 345, 280, 255, 250, 230, 195, 150, 130, 65, 35, 30, 20];
% parasol.density.supportMM.superior = 1:1:22;
% parasol.density.countsMMSq.superior = [1610, 700, 405, 300, 240, 175, 140, 95, 80, 70, 55, 50, 35, 25, 30, 30, 40, nan, nan, nan, nan, nan];
% parasol.density.supportMM.inferior = 1:1:22;
% parasol.density.countsMMSq.inferior = [1820, 1100, 725, 370, 270, 205, 170, 110, 95, 75, 53, 48, 25, 05, nan, nan, nan, nan, nan, nan, nan];
% 
% 


% Parasol proportion of all retinal ganglion cells
% Data from Silviera et al. 1991, Figure 17
parasol.proportion.supportMM.nasal = [0 1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasol.proportion.value.nasal = [0 9.8, 9.64, nan, nan, nan, 7.58, 9.5, 11.74, 13.98, 19.18, 20.78, 15.99, 19.83, 20.32, 24.55, 26.87, 19.93];
parasol.proportion.supportMM.temporal = [0 1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasol.proportion.value.temporal = [0 6.6, 8.44, 9.09, 8.45, 8.85, 10.54, 10.86, 9.51, 9.27, 8.56, 5.69, 2.97, 4.98, nan, nan, nan, nan];
parasol.proportion.supportMM.superior = [0 .96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasol.proportion.value.superior = [0 5.8, 4.71, 4.48, 5.68, 7.91, 8.08, 8.57, 9.38, 9.23, 8.53, 10.44, 11.48, 13.32, 11.9, 9.93];
parasol.proportion.supportMM.inferior = [0 .96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasol.proportion.value.inferior = [0 5.73, 5.58, 7.18, 5.6, 6.96, 8.72, 8.81, 10.01, 8.28, 9.4, 7.74, 9.42, 6.81, 4.68, nan];

% Obtain a spline fit to the parasol densities
for mm = 1:length(p.Results.cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(parasol.proportion.value.(p.Results.cardinalMeridianNames{mm}));
    tmpSupport = parasol.proportion.supportMM.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)';
    parasol.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit( tmpSupport, ...        
        totalRGC.density.fitMMSq.(p.Results.cardinalMeridianNames{mm})(tmpSupport)' .* ...
        parasol.proportion.value.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)' ./100, ...
        'smoothingspline');
end

% Parasol cell body sizes
% Data from Perry et al. 1984, Figure 6C

%% KARA -- The size values below look like microns instead of mm

% Perry, V. H., R. Oehler, and A. Cowey. "Retinal ganglion cells that project to the dorsal
% lateral geniculate nucleus in the macaque monkey." Neuroscience 12.4 (1984): 1101-1123.

parasol.diameter.supportMM = [0.53, 1.07, 1.47, 1.96, 2.5, 3.1, 3.5, 4.0, 4.55, 5.0, 5.53, 6.0, 6.47, 7.0, 7.53, 8.05, 9.0, 9.65, 10.33, 11.4, 12.66, 13.63, 14.21];
parasol.diameter.sizeMM = [0.0113, 0.0126, 0.0138, 0.015, 0.0162, 0.0145, 0.0182, 0.0187, 0.02, 0.0209, 0.0221, 0.0225, 0.0234, 0.0239, 0.0243, 0.0235, 0.0274, 0.026, 0.0242, 0.0238, 0.0227, 0.0224, 0.0201];
parasol.diameter.fitMM = fit(parasol.diameter.supportMM', parasol.diameter.sizeMM','smoothingspline');


%% Bistratified RGCs

% Kara -- you may just have proportion values instead of density for the
% bistratified. Put the numbers here
bistratified.proportion.supportMM.temporal = [];
bistratified.proportion.countsMMSq.temporal = [];

% Bistratified cell body sizes
% Data from Figure 3B:
%   Peterson, Beth B., and Dennis M. Dacey. "Morphology of wide-field
%   bistratified and diffuse human retinal ganglion cells." Visual
%   neuroscience 17.4 (2000): 567-578.

bistratified.diameter.supportMM = [nan];
bistratified.diameter.sizeMM = [0.018];
bistratified.diameter.fitMM = @(x) 0.018;

%% Midget RGCs

% Obtain a spline fit to RGC density values from Curcio & Allen
for mm = 1:length(p.Results.cardinalMeridianAngles)
    midget.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = createMidgetDensityFunc(p.Results.cardinalMeridianAngles(mm));
end


% Midget cell body sizes
% Data from Figure 4B:
%   Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other
%   transparent neurons in the living human retina." Proceedings of the
%   National Academy of Sciences (2017): 201711734.

% Support in the source data is in degrees (of visual field or retina??) 
midget.diameter.supportDeg = [(1.5+3)/2, (3+4.5)/2, (6+7.5)/2, (8+9.5)/2, (12+13.5)/2];
% Are the data from the horizontal meridian? Temporal or nasal retina?
midget.diameter.supportMM = convert_degVisual_to_mmRetina(midget.diameter.supportDeg, 0);
midget.diameter.sizeMM = [0.0115, 0.0113, 0.0114, 0.0118, 0.01315];

% Obtain a spline fit
midget.diameter.fitMM = fit(midget.diameter.supportMM', midget.diameter.sizeMM', 'smoothingspline');

% Figure prep
figure
supportMM = 0:0.01:6;

% Plot counts
subplot(1,2,1)
plot(supportMM, totalRGC.density.fitMMSq.temporal(supportMM) + amacrine.density.fitMMSq.temporal(supportMM))
hold on
plot(supportMM, midget.density.fitMMSq.temporal(supportMM))
plot(supportMM, parasol.density.fitMMSq.temporal(supportMM))
plot(supportMM, amacrine.density.fitMMSq.temporal(supportMM))
plot(supportMM, midget.density.fitMMSq.temporal(supportMM) + parasol.density.fitMMSq.temporal(supportMM) + amacrine.density.fitMMSq.temporal(supportMM),'xr')

sVol = @(d) 4/3*pi*(d./2).^3;

% Plot volume
subplot(1,2,2)
volumeProfile = amacrine.density.fitMMSq.temporal(supportMM) .* sVol(amacrine.diameter.fitMM(supportMM)) + ...
    parasol.density.fitMMSq.temporal(supportMM) .* sVol(parasol.diameter.fitMM(supportMM)./1000) + ...
    midget.density.fitMMSq.temporal(supportMM) .* sVol(midget.diameter.fitMM(supportMM));

plot(supportMM, volumeProfile);

end % rgcThickness function


%% LOCAL FUNCTIONS

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
