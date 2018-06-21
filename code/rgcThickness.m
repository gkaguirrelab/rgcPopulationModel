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
% Curcio & Allen 1990 distinguished displaced amacrine cells from retinal ganglion cells through imagin
% and evalutating their morphology and determined their soma size and densities at eccentricities
% across the human retina. Amacrine cell densities are averages across four meridians.
%
%   Curcio, Christine A., and Kimberly A. Allen. "Topography of ganglion
%   cells in human retina." Journal of comparative Neurology 300.1 (1990):
%   5-25.

% Data from Curcio & Allen 1990, Figure 10
amacrine.density.supportMM.temporal = [0, 0.11, 0.23, 0.4, 0.65, 0.86, 1.45, 2.46, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.7, 18.65, 19.65];
amacrine.density.countsMMSq.temporal = [73, 124, 317, 542, 689, 813, 1052, 1112, 1159, 1187, 1146, 1069, 992, 942, 888, 848, 798, 767, 731, 699, 677, 650, 642, 652, 676];

% Obtain a spline fit to the amacrine densities
for mm = 3:3 %length(p.Results.cardinalMeridianAngles)
    amacrine.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit(amacrine.density.supportMM.(p.Results.cardinalMeridianNames{mm})', amacrine.density.countsMMSq.(p.Results.cardinalMeridianNames{mm})', 'smoothingspline');
end


% Amacrine cell body sizes
% Data from Curcio & Allen 1990, Figure 3
%frequency measurements??
amacrine.diameter.supportMM = [];
amacrine.diameter.sizeMM = [];


%% Parasol RGCs
% Average parasol cell densities across four meridians measured in cells/square mm from six macaque retinas
% Silveira, L. C. L., and V. H. Perry. "The topography of magnocellular projecting ganglion cells (M-ganglion cells)
% in the primate retina." Neuroscience 40.1 (1991): 217-237.

% Data from Silviera et al. 1991, Figures 16
parasol.density.supportMM.temporal = 1:1:22;
parasol.density.countsMMSq.temporal = [2150, 1550, 900, 460, 300, 220, 170, 110, 70, 56, 40, 20, 10, nan, nan, nan, nan, nan, nan, nan, nan, nan];
parasol.density.supportMM.nasal = 1:1:22;
parasol.density.countsMMSq.nasal = [3350, 2560, nan, nan, nan, 520, 485, 415, 400, 380, 345, 280, 255, 250, 230, 195, 150, 130, 65, 35, 30, 20];
parasol.density.supportMM.superior = 1:1:22;
parasol.density.countsMMSq.superior = [1610, 700, 405, 300, 240, 175, 140, 95, 80, 70, 55, 50, 35, 25, 30, 30, 40, nan, nan, nan, nan, nan];
parasol.density.supportMM.inferior = 1:1:22;
parasol.density.countsMMSq.inferior = [1820, 1100, 725, 370, 270, 205, 170, 110, 95, 75, 53, 48, 25, 05, nan, nan, nan, nan, nan, nan, nan];

% Obtain a spline fit to the parasol densities
for mm = 1:length(p.Results.cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(parasol.density.countsMMSq.(p.Results.cardinalMeridianNames{mm}));
    parasol.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit(parasol.density.supportMM.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)', parasol.density.countsMMSq.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)', 'smoothingspline');
end


% Parasol proportion of all retinal ganglion cells
% Data from Silviera et al. 1991, Figure 17
parasol.proportion.supportMM.nasal = [1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasol.proportion.countsMMSq.nasal = [9.8, 9.64, nan, nan, nan, 7.58, 9.5, 11.74, 13.98, 19.18, 20.78, 15.99, 19.83, 20.32, 24.55, 26.87, 19.93];
parasol.proportion.supportMM.temporal = [1.06, 2.07, 3.12, 4.12, 5.09, 6.17, 7.19, 8.22, 9.2, 10.24, 11.26, 12.24, 13.31, 14.32, 15.32, 16.42, 17.39];
parasol.proportion.countsMMSq.temporal = [6.6, 8.44, 9.09, 8.45, 8.85, 10.54, 10.86, 9.51, 9.27, 8.56, 5.69, 2.97, 4.98, nan, nan, nan, nan];
parasol.proportion.supportMM.superior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasol.proportion.countsMMSq.superior = [5.8, 4.71, 4.48, 5.68, 7.91, 8.08, 8.57, 9.38, 9.23, 8.53, 10.44, 11.48, 13.32, 11.9, 9.93];
parasol.proportion.supportMM.inferior = [.96, 1.99, 2.98, 3.93, 4.93, 5.99, 6.87, 7.94, 8.89, 9.88, 10.92, 11.91, 12.9, 13.89, 14.92];
parasol.proportion.countsMMSq.inferior = [5.73, 5.58, 7.18, 5.6, 6.96, 8.72, 8.81, 10.01, 8.28, 9.4, 7.74, 9.42, 6.81, 4.68, nan];

% Parasol cell body sizes
% Data from Perry et al. 1984, Figure 6C

% Perry, V. H., R. Oehler, and A. Cowey. "Retinal ganglion cells that project to the dorsal
% lateral geniculate nucleus in the macaque monkey." Neuroscience 12.4 (1984): 1101-1123.
parasol.diameter.supportMM = [0.53, 1.07, 1.47, 1.96, 2.5, 3.1, 3.5, 4, 4.55, 5, 5.53, 6, 6.47, 7, 7.53, 8.05, 9, 9.65, 10.33, 11.4, 12.66, 13.63, 14.21];
parasol.diameter.sizeMM = [11.27, 12.57, 13.79, 15.02, 16.16, 14.46, 18.23, 18.69, 19.98, 20.9, 22.05, 22.51, 23.35, 23.88, 24.34, 23.5, 27.4, 25.95, 24.19, 23.81, 22.74, 22.35, 20.11];


%% Bistratified RGCs
bistratified.density.supportMM.temporal = [];
bistratified.density.countsMMSq.temporal = [];

% Kara -- you may just have proportion values instead of density for the bistratified. Put the numbers here
bistratified.proportion.supportMM.temporal = [];
bistratified.proportion.countsMMSq.temporal = [];


% Bistratified cell body sizes
% Put comments here about where this info comes from
bistratified.diameter.supportMM = [nan];
bistratified.diameter.sizeMM = [.018];


%% Midget RGCs

% Obtain a spline fit to RGC density values from Curcio & Allen
for mm = 1:length(p.Results.cardinalMeridianAngles)
    midget.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = createMidgetDensityFunc(p.Results.cardinalMeridianAngles(mm));
end


% Midget cell body sizes
% This comes from Liu 2017, Figure 4B
% support in degrees??

% Liu, Zhuolin, et al. "Imaging and quantifying ganglion cells and other transparent neurons
% in the living human retina." Proceedings of the National Academy of Sciences (2017): 201711734.

midget.diameter.supportMM = [1.5-3, 3-4.5, 6-7.5, 8-9.5, 12-13.5];
midget.diameter.sizeMM = [.0115, .0113, .0114, .0118, .01315];


% Calculate thickness

supportMM = 0:0.01:6;
plot(totalRGC.density.fitMMSq.temporal(supportMM) + amacrine.density.fitMMSq.temporal(supportMM))
hold on
plot(midget.density.fitMMSq.temporal(supportMM))
plot(parasol.density.fitMMSq.temporal(supportMM).*100)
plot(amacrine.density.fitMMSq.temporal(supportMM))

plot(midget.density.fitMMSq.temporal(supportMM) + parasol.density.fitMMSq.temporal(supportMM).*100 + amacrine.density.fitMMSq.temporal(supportMM),'xr')

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
