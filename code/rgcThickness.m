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
% Avergae measured parasol cell densities *100 cells/square mm from six macaque retinas
% Silveira, L. C. L., and V. H. Perry. "The topography of magnocellular projecting ganglion cells (M-ganglion cells)
% in the primate retina." Neuroscience 40.1 (1991): 217-237.

% Data from Silviera et al. 1991, Figures 16 and 17
parasol.density.supportMM.temporal = 1:1:22;
parasol.density.countsMMSq.temporal = [21.5, 15.5, 9, 4.6, 3, 2.2, 1.7, 1.1, 0.7, 0.56, 0.4, 0.2, 0.1, nan, nan, nan, nan, nan, nan, nan, nan, nan];
parasol.density.supportMM.nasal = 1:1:22;
parasol.density.countsMMSq.nasal = [33.5, 25.6, nan, nan, nan, 5.2, 4.85, 4.15, 4, 3.8, 3.45, 2.8, 2.55, 2.5, 2.3, 1.95, 1.5, 1.3, 0.65, 0.35, 0.3, 0.2];
parasol.density.supportMM.superior = 1:1:22;
parasol.density.countsMMSq.superior = [16.1, 7, 4.05, 3, 2.4, 1.75, 1.4, 0.95, 0.8, 0.7, 0.55, 0.5, 0.35, 0.25, 0.3, 0.3, 0.4, nan, nan, nan, nan, nan];
parasol.density.supportMM.inferior = 1:1:22;
parasol.density.countsMMSq.inferior = [18.2, 11, 7.25, 3.7, 2.7, 2.05, 1.7, 1.1, 0.95, 0.75, 0.53, 0.48, 0.25, 0.05, nan, nan, nan, nan, nan, nan, nan];

% Obtain a spline fit to the parasol densities
for mm = 1:length(p.Results.cardinalMeridianAngles)
    nonNanSupportIdx = ~isnan(parasol.density.countsMMSq.(p.Results.cardinalMeridianNames{mm}));
    parasol.density.fitMMSq.(p.Results.cardinalMeridianNames{mm}) = ...
        fit(parasol.density.supportMM.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)', parasol.density.countsMMSq.(p.Results.cardinalMeridianNames{mm})(nonNanSupportIdx)', 'smoothingspline');
end


% KARA -- put in the values you digitize from the proportion figure here. You may need to edit the support MM too.
parasol.proportion.supportMM.temporal = 1:1:22;
parasol.proportion.countsMMSq.temporal = [];


% Parasol cell body sizes
% Data from Perry et al. 1984, Figure 6C

% Perry, V. H., R. Oehler, and A. Cowey. "Retinal ganglion cells that project to the dorsal
% lateral geniculate nucleus in the macaque monkey." Neuroscience 12.4 (1984): 1101-1123.
parasol.diameter.supportMM = [];
parasol.diameter.sizeMM = [];


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
